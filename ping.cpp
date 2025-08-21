#include "mui.h"
#include "muiconfig.h"


int main() {
  mui::mpi_split_by_app();
  // Option 1: Declare MUI interface and samplers using specialisms (i.e. 1 = 1 dimensional, d = double)
  printf( "Ping started ....");

  //declare interface
  std::vector<std::unique_ptr<mui::uniface<mui::mui_config>>> mui_ifs;
  // for this case we do not need a spatial sampler since we are 
  // accessing all the individual cells we want
  mui::temporal_sampler_exact<mui::mui_config> temporal_sampler;
  mui::point2d push_point;
  mui::point2d fetch_point;

  std::vector<std::string> ifsName;
  ifsName.emplace_back("ifs1");
  mui_ifs=mui::create_uniface<mui::mui_config>( "ping", ifsName );
	printf("Interface created\n");

  for(int i = 0; i < 100; i++){
    for(int j = 0; j < 100; j++){
        push_point[0] = i;push_point[1] = j;
        mui_ifs[0]->push( "data", push_point, 293.0);
    }
  }

  mui_ifs[0]->commit( 0 );

  printf("Done\n");
  
  // fetch exact values and locations of the points commited this time step
  //std::vector<mui::point<mui::mui_config::REAL, 3>> fetch_locs = mui_ifs[0]->fetch_points<mui::mui_config::REAL>( "dataFromOF", t, temporal_sampler ); // Extract the locations stored in the interface at time=0
  //std::vector<double> fetch_vals = mui_ifs[0]->fetch_values<mui::mui_config::REAL>( "dataFromOF", t, temporal_sampler );
 
  return 0;
}
 