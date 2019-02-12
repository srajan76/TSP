#include <ilcplex/ilocplex.h>
#include <lemon/list_graph.h>
#include <lemon/hao_orlin.h>
#include <lemon/connectivity.h>
#include <optionParser.hpp>
#include "instance.hpp"
#include "model.hpp"
#include "edge.hpp"
#include "hpath.hpp"



int main(int argc, char* argv[]){

    op::OptionParser opt;
    opt.add_option("h", "help", "shows option help" ); 
    opt.add_option("p", "instance_path", "instance_path", "../data/" );
    opt.add_option("f", "file", "name of the instance file", "9602-1-5-0.txt" );
    opt.add_option("n", "targets", "number of targets (other than origin-destination) chosen per TSP", "10" );
    // parse the options and verify that all went well
    bool correct_parsing = opt.parse_options(argc, argv);
    
    if(!correct_parsing) return EXIT_FAILURE;
    if(op::str2bool(opt["h"])) { 
        opt.show_help();
        return 0;
    }

    Instance instance;
    instance.setName(opt["f"]);
    instance.setPath(opt["p"]);
    instance.readData();
     
// The data file from where the coordinates are read has the first row as
// the origin and the last row as the destination (origin == destination)
    
    
    const int n = op::str2int(opt["n"]);
    const int noTargets = instance.getNumTargets();
    const int noofiterations = (noTargets-2)/n;
    std::string path = "../output/";
    std::string filename = path+ "results" + "-n" + std::to_string(n);

    std::ofstream outfile;
    outfile.open(filename);
    outfile<<"Results for TSP size:  " + std::to_string(n) << std::endl;
   for (int k =0; k <noofiterations; ++k) {
    HamiltonianPath hPath(instance);
    hPath.setSource(0);
    hPath.setDestination(noTargets -1);
    int source = hPath.getSource();
    int destination = hPath.getDestination();
    std::vector<int> targetIndexes = {source,destination};
    std::vector<int> satelliteIndexes = {};
        outfile << "target Ids" <<std::endl;
       for ( int j =1+k*n; j < 1+(k+1)*n; ++j){
            targetIndexes.push_back(j);
             outfile <<j<<" ";
        }
    outfile <<std::endl;
    hPath.populatePathData(targetIndexes, satelliteIndexes, source,destination);
    hPath.createEdges();
    hPath.solve();
    outfile <<"OPTIMAL PATH COST FOR TSP"<<std::endl;
    outfile << hPath.getPathCost()<< " " << std::endl;
    outfile <<std::endl;
    outfile <<"OPTIMAL PATH" << std::endl;

    for (int i=0; i<hPath.getPath().size(); ++i) {
        if (i%10 == 0) 
            outfile << std::endl;
        outfile << hPath.getPath()[i] << " ";
    }
    
    outfile << std::endl;

    

   }
    outfile.close();


    return 0;
}



