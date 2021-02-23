
#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector weightedGaussianLocalDensity(NumericVector weights, NumericVector distance, int nrow, double dc) {
  int size = distance.size();
  NumericVector half(size);
  for (int i = 0; i < size; i++) {
    double combOver = distance[i] / dc;
    double negSq = pow(combOver, 2) * -1;
    half[i] = exp(negSq);
  }     
  int ncol = nrow; 
  
  NumericVector result(nrow);
  
  int i = 0;
  for (int col = 0; col < ncol; col++) {
    for (int row = col + 1; row < nrow; row++) {
      double temp = half[i];
      result[row] += temp * (weights[col]);
      result[col] += temp * (weights[row]);
      i++;
    }
  }

  for (int j = 0; j <  nrow; j++){
    result[j] = result[j] + (weights[j] - 1);
  }
  
  return result;
}

// [[Rcpp::export]]
NumericVector weightedNonGaussianLocalDensity(NumericVector weights, int truesize, NumericVector distance, int nrow, double dc) {
  int ncol = nrow;
  NumericVector result(nrow);
  int i = 0;
  for (int col = 0; col < ncol; col++) {
    for (int row = col + 1; row < nrow; row++) {
      if((i % 10000) == 0){
        // if(verbose){
        // Rcout << "index is " << i << " distance under the current index " << distance[i] << std::endl;
        // }
      }
      if(i > distance.size()){
        // Rcout << "Warning: index is larger than the length of the distance vector" << distance[i] << std::endl;
      }
      
      //add in neightbors with weights
      if (distance[i] < dc) {
        result[row] += weights[col];
        result[col] += weights[row];
      } else {
        // do nothing
      }
      i++;
    }
  }

  for (int j = 0; j <  nrow; j++){
    result[j] = result[j] + (weights[j] - 1);
  }
  // if(verbose){
  //  Rcout << "last index is " << i << " length of distance is " << distance.size() << "number of rows is " << nrow << "number of columns is " << ncol << std::endl;
  // }
  return result;
}

// [[Rcpp::export]]
int SumCutOff(NumericVector weights, NumericVector distance, int nrow, double dc) {
  int ncol = nrow;
  int sum = 0;
  int i = 0;
  for (int col = 0; col < ncol; col++) {
    for (int row = col + 1; row < nrow; row++) {
      
      
      //add in neightbors with weights
      if (distance[i] < dc) {
        sum = sum + (weights[col] * weights[row]);
      }
      if (i < (((nrow*(nrow-1))/2)-1)){
        i++;
      }
      
    }
  }
  
  
  // if(verbose){
  //  Rcout << "last index is " << i << " length of distance is " << distance.size() << "number of rows is " << nrow << "number of columns is " << ncol << std::endl;
  // }
  return sum;
}

// [[Rcpp::export]]
NumericVector extractRowColDistfromDistMatrix(NumericVector distance, int num_row, NumericVector row_inds,  NumericVector col_inds) {
  int row_inds_len = row_inds.size();
  int col_inds_len = col_inds.size();
  NumericVector res(row_inds_len * col_inds_len);
  
  // Rcout << "distance is " << distance  << "col_inds index is " << col_inds << std::endl; 
  // Rcout << "row_inds is " << row_inds  << "col_inds index is " << col_inds << std::endl;
  // Rcout << "col_inds_len is " << col_inds_len  << "col_inds_len index is " << col_inds_len  << "res length is " << res.size() << std::endl;
  
  int i = 0;
  int dist_ind; 
  
  for (int row = 0; row < row_inds_len; row++) {
    int row_ind = row_inds[row];
    for (int col = 0; col < col_inds_len; col++) {
      int col_ind = col_inds[col];
      // Rcout << "col is " << col  << "row is " << row << std::endl;
      
      if(row_ind == col_ind){
        res[i] = 0;
      }
      else{
        int row_ind_new; 
        int col_ind_new; 
        
        if(col_ind > row_ind) {
          int row_ind_tmp = row_ind; 
          int col_ind_tmp = col_ind;
          row_ind_new = col_ind_tmp;
          col_ind_new = row_ind_tmp;
        }
        else{
          row_ind_new = row_ind;
          col_ind_new = col_ind;          
        }
        dist_ind = num_row * (col_ind_new - 1) + row_ind_new - 0.5 * (1 + col_ind_new) * col_ind_new - 1;
        // if(row_ind == 3 && col_ind == 2){
        // Rcout << "num_row * (col_ind_new - 1) is " << num_row * (col_ind_new - 1)  << " 1/2 * (1 + col_ind_new) * col_ind_new is " << 0.5 * (1 + 1) * col_ind_new  << " num_row is " << num_row << " dist_ind is " << num_row * (col_ind_new - 1) + row_ind_new - 1/2 * (1 + col_ind_new) * col_ind_new - 1 << std::endl;
        
        // Rcout << "col_ind_new is " << col_ind_new  << " row_ind_new index is " << row_ind_new  << " num_row is " << num_row << " dist_ind is " << dist_ind << " distance under the current index " << distance[dist_ind] << std::endl;
        // }
        res[i] = distance[dist_ind]; 
      }
      i++; 
    }
  }
  
  return res;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

NumericVector all_finite(NumericVector x) {
  return x[x < R_PosInf];
}

// [[Rcpp::export]]
List DensityPeak_smallest_dist_rho_order_coords(NumericVector ordered_rho, NumericVector ordered_coords) {
  int sample_size = ordered_rho.size();
  int dim_num = ordered_coords.size() / sample_size;
  NumericVector smallest_dist(sample_size);
  NumericVector nearest_higher_density_sample(sample_size);
  // Rcout << "sample_size is " << sample_size << "dim_num is " << dim_num << std::endl;
  // Rcout << "ordered_coords is " << ordered_coords << std::endl;
  
  // Rcout << "smallest distances across cells are " << smallest_dist << std::endl;
  
  double current_dist;
  
  for (int cell_ind = 0; cell_ind < sample_size; cell_ind ++) {
    smallest_dist[cell_ind] = R_PosInf;
    nearest_higher_density_sample[cell_ind] = cell_ind;
    
    if(cell_ind == sample_size - 1) { // assign the last distance to the highest density peak cell
      NumericVector all_finite_vals(sample_size - 1);
      all_finite_vals = all_finite(smallest_dist);
      // Rcout << "all_finite_vals are " << all_finite_vals << std::endl;
      NumericVector::iterator maximal =  std::max_element(all_finite_vals.begin(), all_finite_vals.end());
      // Rcout << "maximal index is " << all_finite_vals[maximal - all_finite_vals.begin()] << std::endl;
      smallest_dist[cell_ind] = all_finite_vals[maximal - all_finite_vals.begin()]; 
      nearest_higher_density_sample[cell_ind] = maximal - all_finite_vals.begin();
    }
    
    for (int higher_local_density_cell_ind = cell_ind + 1; higher_local_density_cell_ind < sample_size; higher_local_density_cell_ind ++) {
      // Rcout << "current cell ind is " << cell_ind  << "current cell ind with higher density is " << higher_local_density_cell_ind << std::endl;
      NumericVector source_coord(dim_num);
      NumericVector target_coord(dim_num);
      
      current_dist = 0;
      double tmp;
      for(int dim_num_tmp = 0; dim_num_tmp < dim_num; dim_num_tmp ++) {
        source_coord[dim_num_tmp] = ordered_coords[cell_ind + dim_num_tmp * sample_size];
        target_coord[dim_num_tmp] = ordered_coords[higher_local_density_cell_ind + dim_num_tmp  * sample_size];
        
        tmp = source_coord[dim_num_tmp] - target_coord[dim_num_tmp];
        // tmp = (ordered_coords[cell_ind + 1 + cell_ind * sample_size] - ordered_coords[higher_local_density_cell_ind + 1 + higher_local_density_cell_ind * sample_size]);
        current_dist += tmp * tmp;
      }
      
      current_dist = sqrt(current_dist);
      // Rcout << "current source cell coord is " << source_coord  << "current target cell coord is " << target_coord << std::endl;
      // Rcout << "current_dist is " << current_dist  << std::endl;
      if(smallest_dist[cell_ind] > current_dist) {
        smallest_dist[cell_ind] = current_dist;
        nearest_higher_density_sample[cell_ind] = higher_local_density_cell_ind;
      }
    }
  }
  
  // Rcout << "smallest distances across cells are " << smallest_dist << std::endl;
  // Rcout << "smallest distances across cells are " << nearest_higher_density_sample << std::endl;
  // Rcout << "smallest distances across cells are " << current_dist << std::endl;
  
  // return List::create(Named("smallest_dist") = smallest_dist,
  //                     Named("nearest_higher_density_sample") = nearest_higher_density_sample,
  //                     Named("current_dist") = current_dist);
  return List::create(Named("smallest_dist") = smallest_dist,
                      Named("nearest_higher_density_sample") = nearest_higher_density_sample);
}


/*** R
#extractRowColDistfromDistMatrix(test, attr(test, 'Size'), 1:100, 1:100)
#DensityPeak_smallest_dist_rho_order_coords(rho, master_tsne$Y)
*/

// [[Rcpp::export]]
NumericVector DensityPeak_distancesToPeakC(NumericVector distance, NumericVector rho) {
  int size = rho.size();
  NumericVector peaks(size);
  NumericVector maximum(size);
  
  int i = 0;
  for (int col = 0; col < size; col++) {
    for (int row = col + 1; row < size; row++) {
      double newValue = distance[i];
      double rhoRow = rho[row];
      double rhoCol = rho[col];
      
      if (rhoRow > rhoCol) {
        double peaksCol = peaks[col];
        if (newValue < peaksCol || peaksCol == 0) {
          peaks[col] = newValue;
        }
      } else if (newValue > maximum[col]) {
        maximum[col] = newValue;
      }
      
      if (rhoCol > rhoRow) {
        double peaksRow = peaks[row];
        if (newValue < peaksRow || peaksRow == 0) {
          peaks[row] = newValue;
        }
      } else if (newValue > maximum[row]) {
        maximum[row] = newValue;
      }
      i++;
    }
  }
  
  for (int j = 0; j < size; j++) {
    if (peaks[j] == 0) {
      peaks[j] = maximum[j];
    } else {
      // do nothing, peaks is already min
    }
  } 
  
  return peaks;
}

// [[Rcpp::export]]
NumericVector createWeightedDist(NumericVector distance, int num_row, NumericVector weights) {
  
  int n = num_row;

  NumericVector weightedDist(n*(n-1)/2 + 1);
  
  
  for (int j = 1; j <= n; j++) {
    for (int i = 1; i < j; i++) {
      
      //access previous dissimilarity value between i,j
      int current_index = ((n*(i-1)) - (i*(i-1)/2) + j-i) - 1 ;
      float old_dist = distance[current_index];

      //access the weights of i,j
      int weight_i = weights[i - 1] ;
      int weight_j = weights[j - 1];

      //compute new distance affected by weights
      //equation based off universal gravitation equation  -> NewDist = [(Old Dist)/sqrt(weight_i * weight_j)]
      //float new_dist = old_dist / (sqrt((float)(weight_i * weight_j)));
      //float new_dist = old_dist / (((float)(weight_i * weight_j)));
      //float new_dist = old_dist / (((float)(weight_i + weight_j)/2));
      //float new_dist = old_dist * (((float)(weight_i + weight_j)/2));
      float new_dist = old_dist * (((float)(weight_i * weight_j)));

      //update index in the weighted dist matrix
      weightedDist[current_index] = new_dist;

      


    }
       
  }
  
  return weightedDist;
}

