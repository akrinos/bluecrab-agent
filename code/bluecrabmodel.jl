using Pkg
Pkg.add("Agents")
using Agents

using Pkg
Pkg.add("CSV")
using CSV

global dims_glob = (100,150);
global maxSize = 23; # after this size we no longer grow our crab 
global minTemp = 9; # dormancy temp
global lowTemp = 10.8; # from Brylawski & Miller 2007 - the minimum physiological temp
global lowSalt = 10;
global step = 1;
global nurseryx = 40 # the horizontal limiter on the space for spawn re-entry
global nurseryy = 40 # the vertical limiter on the space for spawn re-entry
println("Initialized.")

global temps_all = CSV.read("model_temps.csv") # the observed temp at each grid cell
global salt_all = CSV.read("model_salt.csv") # the observed salt at each grid cell

println("Read data.")

function getSalt(tempTup,simTime)
    saltstring = salt[convert(Int,tempTup[1]),convert(Int, convert(Int,tempTup[2])+150*simTime)]
    if saltstring == "NA"
        return "NA"
    else
        return parse(Float64, saltstring)
    end
end

function getTemp(tempTup,simTime)
    tempstring = temps[convert(Int,tempTup[1]),convert(Int, convert(Int,tempTup[2])+150*simTime)]
    if tempstring == "NA"
        return "NA"
    else
        return parse(Float64, tempstring)
    end
end

# each agent is a crab that contains its own position, its own ID, its own other physiological params
mutable struct ChesapeakeBayCrab <: AbstractAgent
    id::Int
    pos::NTuple{2, Int}
    temp::Float64 # the temperature where the crab is placed
    salt::Float64 # the salinity where thecrab is placed
    simTime::Int # time since the simulation began
    age::Int
    sex::Int
    size::Float64
    mated::Int
    numSperm::Int
    eggs::Int
    degreeDays::Int
    tOpt::Int
    lifespan::Int
end

# really we have 120 million crabs
function crab_model_2D(;dims = dims_glob, M = 120) 
  global step = 1;
  space = Space(dims, periodic = true)
  age = rand(0:10)
  model = ABM(ChesapeakeBayCrab, space; properties = step, scheduler = random_activation)
  for i in 1:M # add agents in random nodes
      if rand(0:10) > 7
          age = rand(200:1100)
      else
          age = rand(0:200)
      end
      sex = rand(0:1)
      mated = 0
      if rand(1:10) == 1 # fewer need to start out mated. ~10% mated at any given time.
          mated = 1
      end 
      tOpt = rand(20:26)
      degreeDays = rand(0:24)
      simTime = 0
      pos = (rand(1:100), rand(1:150))
      temp = getTemp(pos,simTime)
      while temp == "NA"
          pos = (rand(1:100), rand(1:150))
          temp = getTemp(pos,simTime)
      end
      salt = getSalt(pos,simTime)
      size = rand(0:23)
      numSperm = size * 2
      eggs = round(rand(0:46), digits=0)
      lifespan = rand(200:1100)
      node_to_add = coord2vertex(pos, model)
      add_agent!(node_to_add, model, temp, salt, simTime, age, sex, size, mated, numSperm, eggs, degreeDays, tOpt, lifespan)
  end
  numSpawns = round(0,digits=0)#M/75, digits = 0)
  simTime = 0
  for i in 1:numSpawns
    age = 0
    sex = rand(0:1)
    numSperm = 0
    size = rand(0:5)
    mated = 0
    tOpt = rand(20:26)
    degreeDays = 0
    lifespan = rand(200:1100)
    pos = (rand(nurseryx:100), rand(1:nurseryy)) # spawns can't start just anywhere
    temp = getTemp(pos,simTime)
    numSperm = 0
    eggs = 0
    while temp == "NA"
        pos = (rand(nurseryx:100), rand(1:nurseryy))
        temp = getTemp(pos,simTime)
    end
    salt = getSalt(pos,simTime)
    node_to_add = coord2vertex(pos, model)
    add_agent!(node_to_add, model, temp, salt, simTime, age, sex, size, mated, numSperm, eggs, degreeDays, tOpt, lifespan)
  end
  return model
end

function bay_step!(bay)
    if nagents(bay) > 1
        for node in nodes(bay, by = :random)
            if nagents(bay) == 0
                return
            end
            nc = get_node_contents(node, bay)
            for curr_agent in nc
                crab = id2agent(curr_agent, bay)
                randnumber = rand();
                if crab.age >= crab.lifespan 
                    delete!(bay.agents, crab.id)

                    twoDpos = node
                    splice!(bay.space.agent_positions[twoDpos],
                      findfirst(a->a==crab.id, bay.space.agent_positions[twoDpos]))
                elseif (crab.salt < lowSalt) & (crab.temp < lowTemp) & (randnumber <= (0.25 / 365))
                    delete!(bay.agents, crab.id)
                    twoDpos = node
                    splice!(bay.space.agent_positions[twoDpos],
                      findfirst(a->a==crab.id, bay.space.agent_positions[twoDpos]))
                elseif (crab.temp < lowTemp) & (randnumber <= (0.20 / 365))
                    delete!(bay.agents, crab.id)
                    twoDpos = node
                    splice!(bay.space.agent_positions[twoDpos],
                      findfirst(a->a==crab.id, bay.space.agent_positions[twoDpos]))
                elseif (randnumber <= (0.17 / 365))
                    delete!(bay.agents, crab.id)
                    twoDpos = node
                    splice!(bay.space.agent_positions[twoDpos],
                      findfirst(a->a==crab.id, bay.space.agent_positions[twoDpos]))
                end
            end
        end
    end
end

function agent_step_2d!(agent, model)
    temp = getTemp(agent.pos,agent.simTime)
    salt = getSalt(agent.pos,agent.simTime)
    agent.age += 1
    if (agent.age > 800)
        agent.mated = 2 # don't care about mating anymore
    end
    agent.temp = temp
    agent.salt = salt
    agent.simTime = agent.simTime + 1
    if temp > minTemp
        spawnx = 1 # the horizontal limiter on the space for spawning - any part of lower bay
        spawny = 23 # the vertical limiter on the space for spawning - lower part of bay
        agent_node = coord2vertex(agent.pos, model)
        agent.degreeDays += round(agent.temp - minTemp,digits=0)
        probMolt = 0.669 * log(agent.degreeDays) - 3.71 # from Brylawski & Miller 2007
        growthAmount = real(rand(1044:1598)) / 1000
        if (probMolt > rand()) & ((agent.size + growthAmount) < maxSize)
            agent.size += growthAmount
        end
        if agent.pos[2] > 80
            if agent.sex == 0 # if a female
                agent_node = coord2vertex(agent.pos, model)
                neighboring_nodes = node_neighbors(agent_node, model)
                push!(neighboring_nodes, agent_node) # also consider current node
                available_agents = []
                for curr in neighboring_nodes
                    append!(available_agents, get_node_contents(curr, model))
                end
                done_ = 0
                while length(available_agents) > 0 & done_ == 0
                    tester = rand(available_agents)
                    random_neighbor_agent = id2agent(tester, model)
                    if (random_neighbor_agent.numSperm > 0) & (random_neighbor_agent.sex == 1)
                        agent.mated = 1
                        agent.eggs = random_neighbor_agent.numSperm
                        if (random_neighbor_agent.mated == 1)
                            random_neighbor_agent.mated = 2 # now this sperm donor doesn't care about mating
                            # this will cause there to be a better global distribution of crabs.
                            # and possibly fewer overall crabs.
                        else
                            random_neighbor_agent.mated = 1
                        end
                        done_ = 1
                    end
                
                    deleteat!(available_agents, findfirst(x->x == tester,available_agents))
                end
                
                if (rand(1:1000) == 1) & (agent.mated == 0)
                    agent.eggs = rand(0:23) * 2
                    agent.mated = 1
                end
            else
                agent.numSperm = convert(Int,round(agent.size,digits=0))
            end
        end

        # Randomly decide whether the agent will move position. Each coord can change by up to 1
        newx = agent.pos[1] + rand(-5:5);
        newy = agent.pos[2] + rand(-5:5);

        # Take into account that females that have mated are trying to go southeast (higher x & lower y)
        if (agent.mated == 1) & (agent.sex == 0 | agent.pos[2] > 125)
           newx = agent.pos[1] + rand(-5:5); # + 1 + rand(-5:5);
           newy = agent.pos[2] - 8 + rand(-2:2);
        elseif (agent.pos[2] > 125)
           newy = agent.pos[2] - 4 + rand(-2:2);
        end

        # Take into account that individuals that have not mated are trying to go north & stay west
        if (agent.mated == 0) & (agent.pos[2] < 80) & (agent.pos[1] > 40) & (agent.pos[1] < 65)
            newx = agent.pos[1] + rand(-5:5);
            newy = agent.pos[2] + 8 + rand(-2:2);
        end
        
        # Take into account that individuals that are in branches need to get to the main shaft eventually
        if (agent.mated < 2) & (agent.pos[2] < 80)
            if agent.pos[1] < 40
                newx = agent.pos[1] + rand(-2:4);
            elseif agent.pos[1] > 55
                newx = agent.pos[1] + rand(-4:2); 
            end
        end 
        
        if (agent.pos[1] > 55) & (agent.pos[2] < 40)
            newx = agent.pos[1] - 1 + rand(-4:2);
        end
        
        # distribute crabs that don't care about mating to low-density areas.
        if (agent.mated == 2)
            agent_node = coord2vertex(agent.pos, model)
            currmin = 5
            for cell in node_neighbors(agent_node, model)
                neighboringcrabs = get_node_contents(cell, model)
                if (length(neighboringcrabs) < currmin)
                    (newx,newy) = vertex2coord(cell, model)
                    currmin = length(neighboringcrabs)
                end
            end
        end
        
        # Make sure we don't go off the grid
        if newx > dims_glob[1]
            newx = agent.pos[1]
        elseif newx < 1
            newx = 1
        end

        if newy > dims_glob[2]
            newy = agent.pos[2]
        elseif newy < 1
            newy = 1
        end


        temp = getTemp((newx,newy),agent.simTime)
        numNAs = 0
        while temp == "NA"
            changex = rand(-2:2); 
            changey = rand(-2:2);
            if agent.mated == 0
                changey = rand(-2:5);
            else
                changey = rand(-8:0);
            end
            
            if (agent.mated < 2)
                if agent.pos[1] < 40
                    changex = rand(-1:2);
                elseif agent.pos[1] > 55
                    changex = rand(-2:1); #+ rand(-5:5);
                end
            end
            
            if changex == 0
                changex = rand(-1:1)
            end
            if changey == 0
                changey = rand(-1:1)
            end
            newx = agent.pos[1] + changex 
            newy = agent.pos[2] + changey 
            
            # Make sure we don't go off the grid
            if newx > dims_glob[1]
                newx = dims_glob[1]
            elseif newx < 1
                newx = 1
            end

            if newy > dims_glob[2]
                newy = dims_glob[2]
            elseif newy < 1
                newy = 1
            end
            
            temp = getTemp((newx,newy),agent.simTime)
            numNAs += 1
            
            if (temp == "NA") & (numNAs >= 3)
                newx = agent.pos[1] 
                newy = agent.pos[2]
                temp = agent.temp
            end
            
        end
        

        agent.pos = (newx,newy)
        agent.temp = temp
        agent.salt = getSalt(agent.pos,agent.simTime)
        if (agent.pos[1] > spawnx) & (agent.pos[2] <= spawny) & (agent.mated == 1) & (rand(1:100) == 1)
            numSpawns = round(agent.eggs * 0.3, digits=0) #0.35#rand(20:60); #220);
            simTime = agent.simTime
            for i in 1:numSpawns
                age = 0
                sex = rand(0:1)
                numSperm = 0
                size = rand(0:5)
                mated = 0
                tOpt = rand(20:26)
                degreeDays = 0
                lifespan = rand(200:1100)
                pos = (rand(nurseryx:100), rand(1:nurseryy)) # spawns can't start just anywhere
                temp = getTemp(pos,simTime)
                numSperm = 0
                eggs = 0
                while temp == "NA"
                    pos = (rand(nurseryx:100), rand(1:nurseryy))
                    temp = getTemp(pos,simTime)
                end
                salt = getSalt(pos,simTime)
                node_to_add = coord2vertex(pos, model)
                add_agent!(node_to_add, model, temp, salt, simTime, age, sex, size, mated, numSperm, eggs, degreeDays, tOpt, lifespan) 
            end
            agent.mated = 0
        end
    end
end


using Random
global temps_trunc = CSV.read("model_temps_truncated.csv")[:,1:(150*365)]
global salts_trunc = CSV.read("model_salts_truncated.csv")[:,1:(150*365)]
global temps = hcat(temps_trunc, CSV.read("model_temps_climatechange.csv")[:,1:(150*1460)],makeunique=true) 
# without salt being changed ->
global salt = salt_all[:,1:length(temps[1,:])] 
model = crab_model_2D(;M = 2000)
agent_properties = [:id, :pos, :temp, :simTime, :age, :sex, :mated, :eggs, :size, :numSperm, :degreeDays, :tOpt, :lifespan]
numdays = Int(round(length(temps[1,:]) / 150, digits=0)-1)
data = step!(model, agent_step_2d!, bay_step!, numdays, agent_properties, when =[50,100,150,160,170,180,190,200,300,400,500,600,700,800,900,1000,1100,1400,1700,1800,numdays],replicates=3,parallel=true,step0=true) 

CSV.write(joinpath("output","climatechangetemps_29August2020.csv",data)

global temps = hcat(temps_trunc, CSV.read("model_temps_truncated.csv")[:,1:(150*1460)], makeunique=true) 
model2 = crab_model_2D(;M = 2000)
data2 = step!(model2, agent_step_2d!, bay_step!, numdays, agent_properties, when =[50,100,150,160,170,180,190,200,300,400,500,600,700,800,900,1000,1100,1400,1700,1800,numdays],replicates=3,parallel=true,step0=true) 

Pkg.add("DataFrames")
using DataFrames
CSV.write(joinpath("output","normaltruncatedtemps_29August2020.csv"),data2)


