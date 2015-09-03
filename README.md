
GossipMap is a distributed parallel community detection algorithm to optimize flow-based information-theoretic objective function, called the map equation. 
GossipMap is under GNU General Public License, detailed information is in LICENSE.txt.


## How to compile

GossipMap is implemented in C++ and uses GraphLab PowerGraph for distributed-memory parallelism,
so you have to install GraphLab PowerGraph v2.2 before using GossipMap.
You can find GraphLab PowerGraph from https://github.com/dato-code/PowerGraph.

You can compile GossipMap by following the instruction in the 'Writing Your Own Apps' Section in the GraphLab PowerGraph README.
Below is a modified instruction from the 'Writing Your Own Apps' section for GossipMap application:

    1. Create a sub-directory in the apps/ directory of GraphLab installation, like apps/GossipMap.
    2. Copy GossipMap.cpp and CMakeLists.txt file from the GossipMap directory to apps/GossipMap.
    3. Running 'make' in the apps/ directory should compile GossipMap.  
    4. If GossipMap does not show up, run 'touch apps/CMakeLists.txt' and rerun 'make'


## How to run GossipMap

    1. on a single machine, you could run without using "mpiexec" command.  Also, if you run './GossipMap --help' or './GossipMap' without any arguments, it will show the arguments list with pre-selected values.
        
		* [Usage] >./GossipMap --graph <graph_data> --thresh <threshold> --maxiter <maxIter> --maxspiter <maxSuperIter> --trials <# trials> --mode <1 or 2> --outmode <1 or 2> --ncpus <nCores>
        * [e.g.] >./GossipMap --graph ~/graph_data/web-Stanford.txt --thresh 0.001 --maxiter 10 --maxspiter 3 --trials 1 --mode 1 --outmode 2 --ncpus 8

    2. on multiple machines, you could run GossipMap by using "mpiexec" command.  All of the machines should be installed GraphLab PowerGraph and MPI.
        * [Usage] > mpiexec -f machines /path/to/GossipMap --graph <graph_data> --thresh <threshold> --maxiter <maxIter> --maxspiter <maxSuperIter> --trials <# trials> --mode <1 or 2> --outmode <1 or 2> --ncpus <nCores>

        * 'machines' is a file which contains the hostnames of the machines used for running GossipMap.  
        The command will generate 1 process on each machine represented in 'machines' unless specified.
        If you want to specify the number of MPI processes, you can add '-n <nProcs>' options, such as "mpiexec -n 4 -f machines ..." 
        We recommend use the number of less than or equal to the number of machines for <nProcs> value when you use '-n' option for better performance.

The arguments for GossipMap are following:

  * --help                              
	- Print this help message.
  * --graph arg                         
	- The graph file. Required.
  * --format arg (=snap)                
	- The graph file format. Defaults to (snap). You may use other graphLab readable format.
  * --thresh arg (=0.001)               
	- The threshold for convergence condition. Defaults to (0.001)
  * --tol arg (=1.0000000000000001e-15) 
	- The threshold for pagerank (ergodic state) convergence condition. Defaults to (1E-15)
  * --maxiter arg (=10)                 
	- The maximum of the iteration for finding community. Defaults to (10)
  * --maxspiter arg (=3)                
	- The maximum of the iteration of sp-graph for finding community. Defaults to (3)
  * --trials arg (=1)                   
	- The number of trials for finding community repeatedly. Defaults to (1)
  * --interval arg (=3)                 
	- The time interval for checking whether the received message is valid or not.
  * --mode arg (=1)                     
	- The running mode of finding community: 1 - coreOnce, 2 - coreRepeat. Defaults to (1 = coreOnce).
                                      coreOnce means that GossipMap searches communities with the original-graph once, 
                                      then generate SuperNode graph to search communities with the sp-graph in a SuperStep.
  * --outmode arg (=2)                  
	- The running outerloop mode of finding community: 1 - outerOnce, 2 - outerRepeat. Defaults to (2 = outerRepeat).
                                      'outerOnce' means that GossipMap will run only ONE SuperStep. 
                                      'outerRepeat' will run SuperStep iteratively until GossipMap meet the convergence condition.
  * --prefix arg                        
	- If set, this app will save the community detection result to the given path by the prefix-arg.
                                      If not set, this app will still run and showing log messages but it will not save the community detection result.
  * --ncpus arg (= #cores - 2)          
	- Number of cpus to use per machine. Defaults to (#cores - 2)

There are also some arguments related to GraphLab options.

## Reference

If you would like to add a reference for this application in documents, please put the following bibliography information:

    Seung-Hee Bae and Bill Howe, 
    "GossipMap: A Distributed Community Detection Algorithm for Billion-Edge Directe Graphs,"
    In Proceedings of International Conference on High Performance Computing, Networking, Storage and Analysis (SC'15), 2015 [accepted]


## Contact Information

GossipMap is developed by Seung-Hee Bae and Bill Howe at the University of Washington.
If you want to contact us about GossipMap, you can contact us at:
    
* Seung-Hee Bae: shbae@cs.washington.edu
* Bill Howe: billhowe@cs.washington.edu

Copyright (C) since 2014,  Seung-Hee Bae, Bill Howe, Database Group at the University of Washington
