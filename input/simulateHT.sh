for $i in 1..2
    bin/simulateHighthroughputExperiments/run_simulateHighthroughputExperiments.sh 
        /path/to/runtime 
        seed $i 
        parameterValsPath /path/to/parameterValsPath.{mat|xml} 
        simPath output/sim-$i.mat
end
