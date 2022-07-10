
#include<Rcpp.h>
#include<math.h>

using namespace Rcpp;
// [[Rcpp::export]]
double getDistance(CharacterMatrix flow_path_nodeIDs, NumericVector stream_path_dist,
                   CharacterMatrix total_edgelist) {
  double flow_dist_m = 0.0;
  int flow_row = flow_path_nodeIDs.nrow();

  for(int i=0; i<flow_row; i++) {
    int index = -1;
    for(int j=0; j<total_edgelist.nrow(); j++) {
      if(total_edgelist(j,0)==flow_path_nodeIDs(i,0) &&
	 total_edgelist(j,1)==flow_path_nodeIDs(i,1)) {
	index = j;
	break;
      }
    }
    if(index>0) {
      flow_dist_m += stream_path_dist[index];
    }
  }
  return flow_dist_m;
}
