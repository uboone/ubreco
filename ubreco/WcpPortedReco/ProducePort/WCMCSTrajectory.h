#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TMatrixD.h"
#include "TVector3.h"
#include "TVectorT.h"
#include "TArrayD.h"
#include "TMatrixDEigen.h"

//basic vector functions ==============================================================================================================================================================

std::vector<double> add (std::vector<double> u, std::vector<double> v) {
  for (int i=0,n=u.size();i<n;i++) { u[i] += v[i]; }
  return u;
}

std::vector<double> diff (std::vector<double> u, std::vector<double> v) {
  for (int i=0,n=u.size();i<n;i++) { u[i] -= v[i]; }
  return u;
}

double norm (std::vector<double> v) {
  double norm = 0;
  for (int i=0,n=v.size();i<n;i++) { norm += std::pow(v[i],2); }
  return std::sqrt(norm);
}

double dot (std::vector<double> u, std::vector<double> v) {
  double val = 0;
  for (int i=0,n=u.size();i<n;i++) { val += u[i]*v[i]; }
  return val;
}

double dot_debug (std::vector<double> u, std::vector<double> v) {
  double val = 0;
  for (int i=0,n=u.size();i<n;i++) { val += u[i]*v[i]; 
//  std::cout<<"u1 , u2 , val:  "<<u[i]<<", "<<v[i]<<", "<<val<<std::endl;
  
  }
//  std::cout<<"fixed: "<<(0.942342*0.942054)<<std::endl;
  return val;
}

std::vector<double> cross (std::vector<double> u, std::vector<double> v) {
  return { u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0] };
}

std::vector<double> scale (std::vector<double> v, double a) {
  for (int i=0,n=v.size();i<n;i++) { v[i] *= a; }
  return v;
}

//structs ==============================================================================================================================================================
struct Point;
bool compare_point_score (Point p1, Point p2);
bool compare_edges (std::tuple<double,Point*> t1, std::tuple<double,Point*> t2);

//wrapper structure that lets you sort points with respect to an axis
struct ComparePCAProjection {  
  std::vector<double> axis;
  ComparePCAProjection() {
    (*this).axis = {0,0,0};
  }

  //Helper function for sorting points
  bool compare(std::vector<double> a, std::vector<double> b){ return dot(a,(*this).axis)<dot(b,(*this).axis); }

  //Sort track points in terms of PCA axis
  void sort_points(std::vector<std::vector<double>> &points, std::vector<double> principleAxis){
    (*this).axis = principleAxis;
    //std::sort(points.begin(), points.end(), (*this).compare);
    std::sort(points.begin(), points.end(), [this](std::vector<double> a, std::vector<double> b) { return dot(a,(*this).axis)<dot(b,(*this).axis); });
  }
};

//for use in forming "shortest" path along trajectory points
struct Point {
  int id;
  bool exhausted = false;
  double score;
  std::vector<double> pos, edges_dist;
  Point* prior = NULL;
  std::vector<Point*> edges;

  Point (int id, std::vector<double> pos) {
    (*this).id = id;
    (*this).pos = pos;
    (*this).prior = NULL;
    (*this).score = 1e9;
  }

  //assigns a score to a distance (in cm).  The linear term encourages a short path while the quadratic term demands small distances between segments.
  //These exponents could be changed, particularly if 0.6 cm spacing between points is not used.
  double get_dist_score (Point p) {
    double dist = norm(diff((*this).pos,p.pos));
    return dist + std::pow(dist,2);
  }

  //attempts to update the path to point p if the route through the current point is better than the previous route
  void update_path (int edge_index) {
    Point* edge = (*this).edges[edge_index];
    double new_score = (*this).score + (*this).edges_dist[edge_index];
    if (new_score < (*edge).score) {
      (*edge).score = new_score;
      (*edge).prior = this;
      (*edge).exhausted = false;
    }
  }

  //attempts to update the path to all edges
  void update_paths () {
    for (int i=0,n=(*this).edges.size();i<n;i++) { update_path(i); }
    (*this).exhausted = true;
  }

  void sort_edges () {
    std::vector<std::tuple<double,Point*>> edges_tuples;
    for (int i=0,n=(*this).edges.size();i<n;i++) { edges_tuples.push_back( std::make_tuple((*this).edges_dist[i],(*this).edges[i]) ); }
    std::sort(edges_tuples.begin(), edges_tuples.end(), compare_edges);
    for (int i=0,n=edges_tuples.size();i<n;i++) {
      (*this).edges_dist[i] = std::get<0>(edges_tuples[i]);
      (*this).edges[i]      = std::get<1>(edges_tuples[i]);
    }
  }

  //only add edges along the muon direction (reco end-start) and only keep 4 closest edges
  void add_edge (Point* p, std::vector<double> dir) {
    int nedges_max = 20;
    bool in_dir = dot(diff((*p).pos,(*this).pos), dir) > 0;
    double dist = get_dist_score(*p);
    int nedges = (*this).edges_dist.size();
    if (in_dir && (nedges<nedges_max || dist<(*this).edges_dist.back())) {
      (*this).edges.push_back(p);
      (*this).edges_dist.push_back(dist);
      sort_edges();
      if (nedges>=nedges_max) {
        (*this).edges.pop_back();
        (*this).edges_dist.pop_back();
      }
    }
  }

};

//compares two points based on their score (for use in forming "shortest" path)
struct Comparator {
  bool operator()(const Point* p1, const Point* p2) const { return !(*p1).exhausted && ((*p2).exhausted || ((*p1).score < (*p2).score)); }
};

//compares two point-distance tuples based on their pre-computed distance from a reference point
bool compare_edges (std::tuple<double,Point*> t1, std::tuple<double,Point*> t2) { return std::get<0>(t1) < std::get<0>(t2); }

struct Track {
  std::vector<std::vector<double>> points;
  std::vector<double> weights;
  double total_weight;
  int N; // number of vertices 

  Track () {
    (*this).N            = 0;
    (*this).total_weight = 0;
  }

  void add_point (std::vector<double> point, double weight=1.) {
    (*this).points.push_back(point);
    (*this).weights.push_back(weight);
    (*this).total_weight += weight;
    (*this).N++;
  }

  //removes segment from a track.  Assumes segment and track are sorted the same
  void remove_seg (std::vector<double> first_pos, int size) {
    //find first pos
    int first_index = -1;
    for (int i=0;i<(*this).N;i++) {
      std::vector<double> point_pos = (*this).points[i];
      if (point_pos[0]!=first_pos[0] || point_pos[1]!=first_pos[1] || point_pos[2]!=first_pos[2]) { continue; }
      first_index = i;
      break;
    }

    remove_seg_at(first_index,size);
  }

  //remove a segment of specified size at a specified index
  void remove_seg_at (int first_index, int size) {
    double removed_weight = 0;
    for (int i=0;i<size;i++) { removed_weight += (*this).weights[first_index+i]; } 
    (*this).points.erase( (*this).points.begin()+first_index,  (*this).points.begin() +first_index+size);
    (*this).weights.erase((*this).weights.begin()+first_index, (*this).weights.begin()+first_index+size);
    (*this).total_weight -= removed_weight;
    (*this).N            -= size;

  }

  void clear () {
    (*this).points.clear();
    (*this).weights.clear();
    (*this).total_weight = 0;
    (*this).N = 0;
  }
};

//==============================================================================================================================================================


//Compute the center of mass of a track and set it to the vector com
void setCOM(Track track, std::vector<double> &com) {
  com = {0,0,0};
  for(int i=0; i<track.N; i++){  //compute center of mass
    com[0] += track.points[i][0];
    com[1] += track.points[i][1];
    com[2] += track.points[i][2];
  }
  com = scale(com,1./track.N);
}

//Find the angle separation between two vectors in 3d
double get_angle(std::vector<double> v1, std::vector<double> v2) {
  double l1 = norm(v1);
  double l2 = norm(v2);
  return (acos( dot(v1,v2)/(l1*l2)) );
}

double mat_tot (TMatrixD mat) {
  double tot = 0;
  for (int i=0,m=mat.GetNrows();i<m;i++) {
    for (int j=0,n=mat.GetNcols();j<n;j++) {
      tot += mat(i,j);
    }
  }
  return tot;
}


//PCA Fit ==============================================================================================================================================================

/*
//Computes the RMS displacement for all points in a seg with respect to the (given) pca axis
double getSegRMS(Track seg, std::vector<double> axis, std::vector<double> segCOM){
  double rms = 0;
  std::vector<std::vector<double>> pos;
  TVector3 unitAxis = TVector3(axis[0],axis[1],axis[2]).Unit();
  //Get center of mass coordinates and rotate into coordinates where primary axis is along x, then add to rms
  for(int i=0;i<seg.N;i++){
    TVector3 tempPos = TVector3(seg.points[i][0]-segCOM[0],seg.points[i][1]-segCOM[1],seg.points[i][2]-segCOM[2]);
    tempPos.RotateZ(unitAxis.Phi()*-1);
    tempPos.RotateY(M_PI/2-unitAxis.Theta());
    pos.push_back({tempPos[0],tempPos[1],tempPos[2]});
    rms += pow(pos[i][1],2) + pow(pos[i][2],2);
  }
  rms = sqrt(rms/seg.N);
  return rms;
}
*/

//Helper function that computes the angle between pca fit segments and stores the 3D and 2D angles
void setSegAngles(std::vector<std::vector<double>> priorSegFitVector, std::vector<double> segFitVector,
                  std::vector<double> &angle_vec, std::vector<double> &angleProjX_vec, std::vector<double> &angleProjY_vec){

  std::vector<double> aAxis_prior = priorSegFitVector[0];
  std::vector<double> bAxis_prior = priorSegFitVector[1];
  std::vector<double> cAxis_prior = priorSegFitVector[2];


  std::vector<double> vecy_plane  = cross(aAxis_prior,{1,0,0});
  std::vector<double> vecx_plane  = cross(aAxis_prior,vecy_plane);
  vecx_plane = scale(vecx_plane,1./norm(vecx_plane));
  vecy_plane = scale(vecy_plane,1./norm(vecy_plane));
  std::vector<double> projX = diff(segFitVector, scale(vecy_plane,dot(segFitVector,vecy_plane))) ;
  std::vector<double> projY = diff(segFitVector, scale(vecx_plane,dot(segFitVector,vecx_plane)));
  
  projX = scale(projX,1./norm(projX));
 
  int dirX = 1 - 2*(dot(segFitVector,vecx_plane)<0);
  int dirY = 1 - 2*(dot(segFitVector,vecy_plane)<0);

  angle_vec.push_back(get_angle(segFitVector,aAxis_prior));
  angleProjX_vec.push_back( get_angle(projX, aAxis_prior)*dirX );
  angleProjY_vec.push_back( get_angle(projY, aAxis_prior)*dirY );
}

//Track t is set of 3D points, m is to be computed as mean value of points, the remainder is fit using PCA and stored in a,b,c
std::vector<std::vector<double>> fitPCA(Track track,  std::vector<double> &com, std::vector<double> &evals){
    
  std::vector<std::vector<double>> meanPoints(track.N);
  TArrayD covData(9);
  TMatrixD* covMatrix = new TMatrixD(3,3);
  TVectorD eigenValues;
  evals.resize(3);
	
  com = {0,0,0};
  for(int i=0; i<track.N; i++){  //compute center of mass
    com[0] += track.points[i][0]*track.weights[i];
    com[1] += track.points[i][1]*track.weights[i];
    com[2] += track.points[i][2]*track.weights[i];
  }
  com = scale(com,1./track.total_weight);
  
  for(int i=0; i<track.N; i++){	//create mean value matrix
    meanPoints[i].resize(3);
    meanPoints[i][0] = (track.points[i][0]-com[0]) * track.weights[i];
    meanPoints[i][1] = (track.points[i][1]-com[1]) * track.weights[i];
    meanPoints[i][2] = (track.points[i][2]-com[2]) * track.weights[i];		
  }
	
  //compute and put covariance data into a 1D array
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      for(int k=0; k<track.N; k++) {covData[3*i+j] += meanPoints[k][i] * meanPoints[k][j];}
      covData[3*i+j] = covData[3*i+j] / track.N;
    }
  }
  covMatrix->SetMatrixArray(covData.GetArray());
  
  //find principle axis
  std::vector<std::vector<double>> axes = {{0,0,0},{0,0,0},{0,0,0}};
  const TMatrixD eigenVectors = covMatrix->EigenVectors(eigenValues);

  for (int i=0;i<3;i++){
    axes[0][i] = eigenVectors[i][0];
    axes[1][i] = eigenVectors[i][1];
    axes[2][i] = eigenVectors[i][2];
    evals[i] = eigenValues[i];
  }
  evals = scale(evals, 1./norm(evals));
  delete covMatrix;
  return axes;
}

//helper that takes a track (set of points) and selects a subset of those points to form a track segment, based on the input axis and length.
Track get_seg (Track track, std::vector<double> axis, double seg_length) {
  Track track_seg = Track();
	
  double first_proj = dot(track.points.front(),axis);
  for(int i=0; i<track.N; i++){
    double proj = dot(track.points[i],axis);
    if (proj < (first_proj+seg_length)){
      track_seg.add_point(track.points[i]);
    }
  }

  return track_seg;
}

//Helper Function that fits a segment with all the points in the next seg_length=14 cm starting after the prior lastPointProj
//points is a vector of all the 3D track points
//aAxis is the primary pca axis vector
//seg is a track to be filled with the current selection of points
//segFitVector will store the value of the pca fit for the seg
//segCOM will store the center of mass for the seg
//PointProj variables store info on where to start and end the segments
bool fitSegPCA(Track &track, std::vector<double> aAxis, Track &seg, std::vector<std::vector<double>> &segFitVector, std::vector<double> &segCOM,
               double seg_length, double &currentFirstPointProj, double &currentLastPointProj){
  //setup
  seg.points.clear();
  seg.N = 0;
	
  //track_buffer stores only a short amount of the track to speed up computation when sorting
  Track track_buffer = get_seg(track, aAxis, seg_length*2);

  //update currentFirstPointProj
  double lastPointProj = dot(track_buffer.points.back(),aAxis);
  currentFirstPointProj = lastPointProj;
  for(int i=0; i<int(track_buffer.points.size()); i++){
    double proj = dot(track_buffer.points[i],aAxis);
    if((proj>currentLastPointProj) && (proj<currentFirstPointProj)){ currentFirstPointProj = proj; }
  }

  //fill seg with points
  seg = get_seg(track_buffer, aAxis, seg_length);
  if (seg.N<=1) { return false; }

  //run PCA on the segment and orient along aAxis
  std::vector<double> eigenvalues;
  segFitVector = fitPCA(seg, segCOM, eigenvalues);
  if(dot(segFitVector[0],aAxis)<0) {
    segFitVector[0][0] *= -1;
    segFitVector[0][1] *= -1;
    segFitVector[0][2] *= -1;	
  }

  //Sort points along new fit axis
  ComparePCAProjection comparator = ComparePCAProjection();
  comparator.sort_points(track_buffer.points, segFitVector[0]);

  //re-determine what points belong in the segment based on new u-vector
  seg = get_seg(track_buffer, segFitVector[0], seg_length);
  if (seg.N<=1) { return false; }

  //refit and orient along aAxis
  segFitVector = fitPCA(seg, segCOM, eigenvalues);
  if(dot(segFitVector[0],aAxis)<0) {
    segFitVector[0][0] *= -1;
    segFitVector[0][1] *= -1;
    segFitVector[0][2] *= -1;	
  }

  //update current-last-point-proj
  comparator.sort_points(seg.points, segFitVector[0]);
  currentLastPointProj = dot(seg.points.back(),aAxis);

  //remove points from original track
  comparator.sort_points(seg.points, aAxis);
  track.remove_seg(seg.points.front(),seg.N);
  return true;
}

//Helper function to track the largest and smallest track positions along the b and c axis for setting up graph bounds
void setTrackBounds(std::vector<double> &bounds, std::vector<double> bProj, std::vector<double> cProj, std::vector<double> bAxis, std::vector<double> cAxis){
  for(int i=0;i<int(bProj.size());i++){
    if (bProj[i]<bounds[0] || bounds[0]==0) {bounds[0] = bProj[i];}
    if (bProj[i]>bounds[1] || bounds[1]==0) {bounds[1] = bProj[i];}
    if (cProj[i]<bounds[2] || bounds[2]==0) {bounds[2] = cProj[i];}
    if (cProj[i]>bounds[3] || bounds[3]==0) {bounds[3] = cProj[i];}
  }
}



//helper function that trims an input list of trajectory points to only include those along the shortest path from muon start to end vertices
std::tuple<bool,std::vector<std::vector<double>>> trim_trajectory(double npoints_traj, std::vector<std::vector<double>> trajectory_points_initial, std::vector<double> vtx_muon_start_reco, std::vector<double> vtx_muon_end_reco) {

  //trim trajectory points to get shortest path along muon track
  //find starting and ending trajectory points based on proximinty to reco muon start and end 
  std::vector<double> muon_dir_reco = diff(vtx_muon_end_reco,vtx_muon_start_reco);
  int startpoint_index = -1;
  int endpoint_index = -1;
  double startpoint_dist = 1e9;
  double endpoint_dist = 1e9;
  for (int i=0;i<npoints_traj;i++) {
    std::vector<double> pos = trajectory_points_initial[i];
    double dist_start = std::sqrt( norm(diff(pos,vtx_muon_start_reco)) );
    double dist_end   = std::sqrt( norm(diff(pos,vtx_muon_end_reco)) );
    if (dist_start < startpoint_dist) {
      startpoint_index = i;
      startpoint_dist = dist_start;
    }
    if (dist_end < endpoint_dist) {
      endpoint_index = i;
      endpoint_dist = dist_end;
    }
  }
  //create a vector of all trajectory points, with the 0th Point as the first in the vector. Give the starting point a graph score of 0.
  std::vector<Point*> traj_points;
  Point* new_point0 = new Point(startpoint_index+1, { trajectory_points_initial[startpoint_index][0], trajectory_points_initial[startpoint_index][1], trajectory_points_initial[startpoint_index][2] });
  traj_points.push_back(new_point0);
  (*traj_points[0]).score = 0;
  for (int i=0;i<npoints_traj;i++) {
    Point* new_point = new Point(i+1, { trajectory_points_initial[i][0], trajectory_points_initial[i][1], trajectory_points_initial[i][2] });
    if (i != startpoint_index) { traj_points.push_back(new_point); }
  }
  //form graph of all Points by taking each Point i and calling add_edge on each other point j
  for (int i=0;i<npoints_traj;i++) { for (int j=0;j<npoints_traj;j++) { (*traj_points[i]).add_edge(traj_points[j],muon_dir_reco); } }
  //call update_paths on 0th Point in the vector
  //sort points by their score
  //check whether the 0th Point in the vector is the end point.  If not, loop again
  bool bad_path = false;
  while (true) {
    (*traj_points.front()).update_paths();
    std::sort(traj_points.begin(), traj_points.end(), Comparator());
    if ((*traj_points.front()).id == (endpoint_index+1)) { break; }
    if ((*traj_points.front()).exhausted) {
      std::cout << "unable to traverse particle path" << std::endl;
      bad_path = true;
      break;
    }
  }
  //remove all trajectory Points that are not in the path to the end point
  std::vector<std::vector<double>> trajectory_points_final;
  Point* path_pointer = traj_points[0];
  while(true && !bad_path) {
    trajectory_points_final.push_back((*path_pointer).pos);
    path_pointer = (*path_pointer).prior;
    if (path_pointer == NULL) { break; }
  }

  //add muon start and end vertex to list of trajectory points
  trajectory_points_final.insert(trajectory_points_final.begin(), vtx_muon_end_reco);
  trajectory_points_final.insert(trajectory_points_final.end(),   vtx_muon_start_reco);

  //clean up
  for (int i=0;i<npoints_traj;i++) { delete traj_points[i]; }

  return std::make_tuple(bad_path, trajectory_points_final);
}


//take in point cloud, muon start, and specified segment length
//Forms segments and fits each
std::tuple< std::vector<Track>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double> > form_segs (std::vector<std::vector<double>> vec_points, std::vector<double> muon_start, std::vector<double> muon_end, double seg_length) {
  std::vector<Track> seg_vec;            //vector of each segment of trach
  std::vector<double> distance_vec;      //distance from muon start to middle of each segment
  std::vector<double> angle_vec;         //angle between each segment
  std::vector<double> displacement_vec;  //rms displacement of points within each segment

  //fill track
  Track track;
  int npoints = vec_points.size();
  for(int i=0; i<npoints; i++){ track.add_point(vec_points[i]); }

  //get PCA axes for track
  std::vector<double> com         = {0,0,0};
  std::vector<double> eigenvalues = {0,0,0};
  std::vector<std::vector<double>> axes = fitPCA(track, com, eigenvalues);
  std::vector<double> aAxis = axes[0];
  std::vector<double> bAxis = axes[1];
  std::vector<double> cAxis = axes[2];
  //flip axis to be along muon direction
  std::vector<double> vec_muon = { muon_end[0]-muon_start[0], muon_end[1]-muon_start[1], muon_end[2]-muon_start[2] };
  if(dot(vec_muon,aAxis) < 0){ aAxis = scale(aAxis,-1.); }
  ComparePCAProjection comparator = ComparePCAProjection();
  comparator.sort_points(track.points, aAxis);

  //fill temporary track with remaining points, this gets reduced as they are added to segments
  Track track_remainder = get_seg(track, aAxis, 1e6);

  std::vector<std::vector<std::vector<double>>> segPCAFit_vec;  //pca fit vectors A,B,C for each 14 cm segment
  std::vector<std::vector<double>> segs_aAxis_vec;  //center of mass position for each segment
  std::vector<std::vector<double>> segCOM_vec;  //center of mass position for each segment
  std::vector<double> angleProjB_vec, angleProjC_vec;

  double currentFirstPointProj = dot(track.points.front(),aAxis);
  double currentLastPointProj  = dot(track.points.front(),aAxis);
  double end_proj              = dot(track.points.back(), aAxis);

  //subdivide into segments of specified length
  //use flexible length to handle gaps
  int iseg = 0;
  while((currentLastPointProj+0.5*seg_length) < end_proj && track_remainder.N>=10){

    segPCAFit_vec.push_back({{0,0,0}, {0,0,0}, {0,0,0}});
    segCOM_vec.push_back(   {0,0,0});
    seg_vec.push_back(Track());

    //Form and fit 14 cm seg with pcafit, update currentFirstPointProj, currentLastPointProj
    bool can_fit = fitSegPCA(track_remainder, aAxis, seg_vec.back(), segPCAFit_vec.back(), segCOM_vec.back(), seg_length, currentFirstPointProj, currentLastPointProj);
    if (!can_fit) { break; } 
    segs_aAxis_vec.push_back(segPCAFit_vec.back()[0]);
    distance_vec.push_back((currentFirstPointProj+currentLastPointProj)/2 - dot(muon_start,aAxis));

    //get 3D angle between tracks and 2D projection angles
    if (iseg==0) { angle_vec.push_back(-1); angleProjB_vec.push_back(-1); angleProjC_vec.push_back(-1); }
    else         { setSegAngles(segPCAFit_vec[iseg-1], segs_aAxis_vec.back(), angle_vec, angleProjB_vec, angleProjC_vec); }

    iseg++;
  }

  return { seg_vec, axes, segs_aAxis_vec, segCOM_vec, distance_vec, angle_vec, angleProjB_vec, angleProjC_vec };
}


/*
//Point cloud methods (deprecated) ==============================================================================================================================================================

//Takes in a vector of points as well as muon start,end reco coordinates and returns a subset of those points that belong to the muon
//aimed at rejecting non-muon particles (proton, etc) and crossing tracks
std::vector<std::vector<std::vector<double>>> select_points_muon_simple (std::vector<std::vector<double>> vec_points_all, std::vector<double> muon_start, std::vector<double> muon_end, std::vector<std::vector<double>> input_simple_path ) {

  double threshold_angle      = 30;       //degrees, reject points that are not along the muon direction
  double threshold_start      = 5;        //cm, reject points near start of muon track
  double threshold_end        = 2;        //cm, reject points near end of muon track
  double threshold_perp_const = 40;       //reject points not along line from muon start to end
  double threshold_perp_scale = 2.0/1000; //increase tolerance by 0.5cm for every meter of track
  double segment_length       = 2.0;      //length of each segment subdivision of points along muon vector
  double threshold_path       = 2.0;      // 2.8  threshold for selecting points near shortest path through graph

  //reject points that are not along a line from muon start to end
  std::vector<std::vector<double>> vec_points_candidate;
  int npoints_all = vec_points_all.size();
  std::vector<double> vec_muon = { muon_end[0]-muon_start[0], muon_end[1]-muon_start[1], muon_end[2]-muon_start[2] };
  double vec_muon_length = norm(vec_muon);
  vec_muon = scale(vec_muon, 1./vec_muon_length);
  double mu_start_dot = dot(muon_start,vec_muon);
  double mu_end_dot   = dot(muon_end,vec_muon);

  for (int i=0;i<npoints_all;i++) {
    std::vector<double> point_pos = vec_points_all[i];
    std::vector<double> point_pos_adjusted = diff(point_pos, muon_start);
    double point_mu_angle = (180./M_PI)*std::acos(dot(point_pos_adjusted,vec_muon)/norm(point_pos_adjusted));
    double point_mu_dot = dot(point_pos,vec_muon);
    //require points between muon start and end and along muon direction
    if ( point_mu_angle>threshold_angle || point_mu_dot<mu_start_dot+threshold_start || point_mu_dot>mu_end_dot-threshold_end ) { continue; }

    //require points that are along a wide linear path from muon start to end
    std::vector<double> point_mu_parallel = scale(vec_muon, dot(point_pos_adjusted,vec_muon));
    std::vector<double> point_mu_perp     = diff(point_pos_adjusted,point_mu_parallel);
    if (norm(point_mu_perp) > threshold_perp_const+threshold_perp_scale*vec_muon_length) { continue; }

    vec_points_candidate.push_back(vec_points_all[i]);
  }

  int npoints_candidate = vec_points_candidate.size();
  if (npoints_candidate==0) { return {{},{},{}};}



  //group points into segments along projection from start to end
  sort_points(vec_points_candidate, vec_muon);
  std::vector<std::vector<std::vector<double>>> vec_points_segmented = {{}};
  double old_point_proj = dot(vec_points_candidate[0],vec_muon);
  int seg_index = 0;
  for (int i=0;i<npoints_candidate;i++) {
    double proj = dot(vec_points_candidate[i],vec_muon);
    if (proj>old_point_proj+segment_length) {
      old_point_proj = proj;
      vec_points_segmented.push_back({});
      seg_index++;
    }
    vec_points_segmented[seg_index].push_back(vec_points_candidate[i]);
  }

  //find point with lowest offset and add it to simple path vector
  std::vector<double> prior_seg_pos = muon_start;
  std::vector<double> vec_point_end = vec_muon;
  std::vector<std::vector<double>> vec_points_path;
  if (!input_simple_path.empty()) { vec_points_path = input_simple_path; }
  else {
    for (int i=0,nsegs=vec_points_segmented.size();i<nsegs;i++) {
      double min_offset = 1e9;
      int min_offset_index = -1;
      for (int j=0,m=vec_points_segmented[i].size();j<m;j++) {
        std::vector<double> point_pos_adjusted = diff(vec_points_segmented[i][j], prior_seg_pos);
        std::vector<double> point_pos_parallel = scale(vec_point_end, dot(point_pos_adjusted,vec_point_end));
        std::vector<double> point_pos_perp     = diff(point_pos_adjusted, point_pos_parallel);
        double offset = norm(point_pos_perp);
        if (offset<min_offset) {
          min_offset = offset;
          min_offset_index = j;
        }
      }

      std::vector<double> min_offset_point = vec_points_segmented[i][min_offset_index];
      vec_points_path.push_back(min_offset_point);
      prior_seg_pos = min_offset_point;
      vec_point_end = diff(muon_end,min_offset_point);
      vec_point_end = scale(vec_point_end, 1./norm(vec_point_end));
    }
  }
    
  //reject points not along simple path
  std::vector<bool> vec_near_path(npoints_candidate, false);
  for (int i=0,n=vec_points_path.size();i<n;i++) {
    for (int j=0,m=vec_points_candidate.size();j<m;j++) {
      if (norm(diff(vec_points_path[i], vec_points_candidate[j])) < threshold_path) { vec_near_path[j] = true; }
    }
  }
  std::vector<std::vector<double>> vec_points_selected;
  for (int i=0;i<npoints_candidate;i++) { if (vec_near_path[i]) { vec_points_selected.push_back(vec_points_candidate[i]); } }
  int npoints_selected = vec_points_selected.size();

  return { vec_points_all, vec_points_candidate, vec_points_selected };
}

//Helper function that computes COM, cov, fits tangent vector, and computes updated COM
std::vector<std::vector<double>> compute_averaged_COM (std::vector<std::vector<std::vector<double>>> v_segs_points, std::vector<std::vector<double>> &v_segs_tangent, std::vector<double> aAxis, std::vector<double> first_point, double seg_length, int local_range_index) {
  
  int nsegs = v_segs_points.size();
  std::vector<std::vector<double>> v_seg_COM(        nsegs,{0,0,0}); //center of mass position for each segment
  std::vector<std::vector<double>> v_seg_COM_updated(nsegs,{0,0,0}); //update center of mass position for each segment based on nearby segments
  std::vector<TMatrixD> vm_seg_covCOM;

  //compute center of mass and uncertainty for each segment
  for (int i=0;i<nsegs;i++) {
    std::vector<double> seg_COM = {0,0,0};
    double npoints_seg = v_segs_points[i].size();
    TArrayD covData(9);
    TMatrixD m_covCOM(3,3);
    std::vector<std::vector<double>> v_seg_mean_points;
    //center of mass
    if (npoints_seg>0) {
      for (int j=0;j<npoints_seg;j++) { seg_COM = add(seg_COM,v_segs_points[i][j]); }
      seg_COM = scale(seg_COM, 1./npoints_seg);
      for (int j=0;j<npoints_seg;j++) { v_seg_mean_points.push_back(diff(v_segs_points[i][j],seg_COM)); }
    } else {
      if (i>0) { seg_COM = add(v_seg_COM[i-1], scale(aAxis,seg_length)); }
      else     { seg_COM = first_point; }
    }
    //covariance
    for (int j=0;j<3;j++) {
      for (int k=0;k<3;k++) {
        for (int ii=0;ii<npoints_seg;ii++) { covData[3*j+k] += v_seg_mean_points[ii][j]*v_seg_mean_points[ii][k]; }
          covData[3*j+k] /= npoints_seg;          
        }
    }
    if (npoints_seg<3 || mat_tot(m_covCOM)==0) { for (int j=0;j<9;j++) { covData[j] = 30; } }
    v_seg_COM[i]      = seg_COM;
    m_covCOM.SetMatrixArray(covData.GetArray());
    vm_seg_covCOM.push_back(m_covCOM);
  }

  //fit tangent line at each point using nearby points
  //Then, update COM for target point by taking average of points projected onto plane perpendicular to tangent line
  for (int i=0;i<nsegs;i++) {
    int min_index = std::clamp(i-local_range_index,0,nsegs-1);
    int max_index = std::clamp(i+local_range_index,0,nsegs-1);
    Track t_local_points;
    for (int j=min_index;j<=max_index;j++) { t_local_points.add_point( v_seg_COM[j],1./sqrt(mat_tot(vm_seg_covCOM[j])) ); }
    std::vector<double> seg_COM     = {0,0,0};
    std::vector<double> eigenvalues = {0,0,0};
    //get PCA fit vectors A,B,C.  A is tangent line
    std::vector<std::vector<double>> fitvecs = fitPCA(t_local_points,seg_COM,eigenvalues);
    v_segs_tangent[i] = fitvecs[0];
    TMatrixD m_fitvec_B(3,1);
    TMatrixD m_fitvec_C(3,1);
    for (int j=0;j<3;j++) {
      m_fitvec_B[j][0] = fitvecs[1][j];
      m_fitvec_C[j][0] = fitvecs[2][j];
    }

    TMatrixD m_fitvec_B_t(1,3);
    TMatrixD m_fitvec_C_t(1,3);
    m_fitvec_B_t.Transpose(m_fitvec_B);
    m_fitvec_C_t.Transpose(m_fitvec_C);
    //project COM onto A,B,C axes.  Update B,C with contributions from nearby points
    std::vector<double> seg_COM_projA = scale(fitvecs[0], dot(fitvecs[0],v_seg_COM[i]));
    std::vector<double> seg_COM_projB = {0,0,0};
    std::vector<double> seg_COM_projC = {0,0,0};
    double sigmaB_total = 0;
    double sigmaC_total = 0;
    for (int j=min_index;j<=max_index;j++) {
      double sigmaB = sqrt((m_fitvec_B_t*vm_seg_covCOM[j]*m_fitvec_B)[0][0]);
      double sigmaC = sqrt((m_fitvec_C_t*vm_seg_covCOM[j]*m_fitvec_C)[0][0]);
      sigmaB_total += 1./sigmaB;
      sigmaC_total += 1./sigmaC;
      seg_COM_projB = add(seg_COM_projB, scale( fitvecs[1], dot(fitvecs[1],v_seg_COM[j])/sigmaB));  
      seg_COM_projC = add(seg_COM_projC, scale( fitvecs[2], dot(fitvecs[2],v_seg_COM[j])/sigmaC ));
    }
    seg_COM_projB = scale(seg_COM_projB, 1./sigmaB_total);
    seg_COM_projC = scale(seg_COM_projC, 1./sigmaC_total);

    //update COM
    v_seg_COM_updated[i] = add(v_seg_COM_updated[i],seg_COM_projA);
    v_seg_COM_updated[i] = add(v_seg_COM_updated[i],seg_COM_projB);
    v_seg_COM_updated[i] = add(v_seg_COM_updated[i],seg_COM_projC);
  }
  return v_seg_COM_updated;

}

//Helper function that trims points from a segmented point cloud that are not close to a set of select points (one per seg required) by radius
std::vector<std::vector<std::vector<double>>> trim_seg_points_radius (std::vector<std::vector<std::vector<double>>> seg_points, std::vector<std::vector<double>> select_points, double threshold) {

  //label points as kept or rejected
  int nsegs = seg_points.size();
  std::vector<std::vector<bool>> v_seg_points_near_path;
  for (int i=0;i<nsegs;i++) {
    int npoints = seg_points[i].size();
    std::vector<bool> v_points_near_path(npoints, false);
    int min_index = std::clamp(i-1,0,nsegs-1);
    int max_index = std::clamp(i+1,0,nsegs-1);
    for (int j=min_index;j<=max_index;j++) {
      for (int k=0;k<npoints;k++) {
        if (norm(diff(seg_points[i][k],select_points[j]))<threshold) { v_points_near_path[k] = true; }
      }
    }
    v_seg_points_near_path.push_back(v_points_near_path);
  }

  //create output
  std::vector<std::vector<std::vector<double>>> v_seg_points_selected;
  for (int i=0;i<nsegs;i++) {
    v_seg_points_selected.push_back({});
    int npoints = seg_points[i].size();
    for (int k=0;k<npoints;k++) {
      if (v_seg_points_near_path[i][k]) { v_seg_points_selected[i].push_back(seg_points[i][k]); }
    }
  }
  return v_seg_points_selected;

}

//helper function that takes in two points and vectors that form a plane, and computes the distance between the points projected onto the plane
double get_plane_offset (std::vector<double> p, std::vector<double> q, std::vector<double> perp_a, std::vector<double> perp_b) {
  std::vector<double> p_proj = {0., dot(perp_a,p), dot(perp_b,p)};
  std::vector<double> q_proj = {0., dot(perp_a,q), dot(perp_b,q)};
  return norm(diff(p_proj,q_proj));
}

//Helper function that computes the perpendicular offset w.r.t. a centerpoint and given direction that contians 90% of points in input pointcloud
double get_offset_90p (std::vector<std::vector<double>> v_points, std::vector<double> dir, std::vector<double> reference_dir, std::vector<double> center, std::vector<double> &center_perp, std::vector<double> &perpU, std::vector<double> &perpV, bool compute_perp) {
  if (v_points.size()==0) { return 0; }
  if (compute_perp || perpU.size()==0) {
    perpU = cross(dir,reference_dir);
    perpV = cross(dir,perpU);
    perpU = scale(perpU,1./norm(perpU));
    perpV = scale(perpV,1./norm(perpV));
  }

  center_perp = { 0, dot(center,perpU) , dot(center,perpV) };
  std::vector<double> v_offset;
  for (int i=0,n=v_points.size();i<n;i++) {
    std::vector<double> point_perp = { 0, dot(v_points[i],perpU) , dot(v_points[i],perpV) };
    v_offset.push_back(norm(diff(point_perp,center_perp)));
  }
  std::sort(v_offset.begin(), v_offset.end());
  int offset_90p_index = (int)(v_offset.size()*9./10);
  return v_offset[offset_90p_index];
}

//Helper function that trims points from a segmented point cloud that are not close to a set of select points (one per seg required) by perpendicular offset
//Averages the measured max offset across local segments
std::vector<std::vector<std::vector<double>>> trim_seg_points_offset (std::vector<std::vector<std::vector<double>>> seg_points, std::vector<std::vector<double>> select_points, std::vector<std::vector<double>> v_segs_tangent, std::vector<double> &v_segs_offset_90p_ave, int local_range_index, double threshold_ratio) {

  //for each seg, get perpendicular offset (using tangent line) between points and the COM point for that segment, find value that contains 90% of points
  int nsegs = seg_points.size();
  std::vector<double> v_segs_offset_90p;
  std::vector<std::vector<double>> v_segs_perp_a, v_segs_perp_b;
  for (int i=0;i<nsegs;i++) {
    std::vector<double> select_point_proj, perp_a, perp_b;
    double offset_90p = get_offset_90p(seg_points[i], v_segs_tangent[i], {1.,0.,0.}, select_points[i], select_point_proj, perp_a, perp_b, true);
    v_segs_offset_90p.push_back(offset_90p);
    v_segs_perp_a.push_back(perp_a);
    v_segs_perp_b.push_back(perp_b);
  }

  //Update 90% offset values by averaging across nearby points
  for (int i=0;i<nsegs;i++) {
    v_segs_offset_90p_ave[i] = 0;
    int min_index = std::clamp(i-local_range_index,0,nsegs-1);
    int max_index = std::clamp(i+local_range_index,0,nsegs-1);
    int nsegs_nonzero = 0;
    for (int j=min_index;j<=max_index;j++) {
      if (v_segs_offset_90p[j]>0) {
        nsegs_nonzero++;
        v_segs_offset_90p_ave[i] += v_segs_offset_90p[j];
      }
    }
    v_segs_offset_90p_ave[i] /= std::max(1,nsegs_nonzero);
  }

  //label points as kept or rejected based on perpendicular offset from selected (COM) point
  std::vector<std::vector<bool>> v_seg_points_near_path;
  for (int i=0;i<nsegs;i++) {
    int npoints = seg_points[i].size();
    std::vector<bool> v_points_near_path(npoints, false);
    for (int k=0;k<npoints;k++) {
      double offset = get_plane_offset(seg_points[i][k], select_points[i], v_segs_perp_a[i], v_segs_perp_b[i]);
      if (offset<v_segs_offset_90p_ave[i]*threshold_ratio + 0.1) { v_points_near_path[k] = true; } //added 0.1 cm allowance
    }
    v_seg_points_near_path.push_back(v_points_near_path);
  }

  //create output
  std::vector<std::vector<std::vector<double>>> v_seg_points_selected;
  for (int i=0;i<nsegs;i++) {
    v_seg_points_selected.push_back({});
    int npoints = seg_points[i].size();
    for (int k=0;k<npoints;k++) {
      if (v_seg_points_near_path[i][k]) { v_seg_points_selected[i].push_back(seg_points[i][k]); }
    }
  }
  return v_seg_points_selected;
}

//Helper function that takes a prong COM and projectes it onto the center of the main track
std::vector<double> generate_track_point (std::vector<double> prong_COM, std::vector<std::vector<double>> v_segs_COM, std::vector<std::vector<double>> v_segs_tangent, std::vector<double> aAxis) {
  double prong_COM_proj = dot(aAxis,prong_COM);
  std::vector<double> v_seg_dist;
  for (int i=0,n=v_segs_COM.size();i<n;i++) {
    double seg_COM_proj = dot(aAxis,v_segs_COM[i]);
    v_seg_dist.push_back(abs(prong_COM_proj-seg_COM_proj));
  }
  int seg_index = std::distance(v_seg_dist.begin(),std::min_element(v_seg_dist.begin(), v_seg_dist.end()));
  std::vector<double> seg_COM = v_segs_COM[seg_index];
  std::vector<double> seg_tangent = v_segs_tangent[seg_index];
  std::vector<double> seg_COM_perp = diff( seg_COM, scale(seg_tangent,dot(seg_tangent,seg_COM)) );
  prong_COM_proj = dot(seg_tangent,prong_COM);
  return add( seg_COM_perp, scale(seg_tangent,prong_COM_proj) );
}

//Helper function that identifies points that are part of a prong (crossing track, delta ray) and removes them
//Takes in track of previously rejected points, vector of segs of points, vector of seg COMs, vector of seg tangents, vector of seg 90% offsets, aAxis of overall point cloud, length of each seg (14cm), and projection of first and last points along aAxis 
//returns the prong COMs and dirs
std::vector<std::vector<std::vector<double>>> remove_prongs (Track &track_reject, std::vector<std::vector<std::vector<double>>> &v_segs_points_preliminary, std::vector<std::vector<double>> averaged_COM_preliminary, std::vector<std::vector<double>> v_segs_tangent, std::vector<double> v_segs_offset_90p, std::vector<double> aAxis, double seg_length, double firstPointProj, double lastPointProj, std::vector<double> bAxis, std::vector<double> cAxis) {

  int npoints_seg_reject_threshold   = 16;
  double seg_reject_length           = 6.0;  //determines ROI for finding a prong
  double prong_width_threshold       = 3.0;  //allowed width for points to be considered part of prong
  double prong_width_ratio_threshold = 1.2;  //only label prong points within a ratio of the 90-percentile effective prong radius
  double prong_angle_threshold       = 12.0*M_PI/180;  //minimum angle required to be a prong
  double prong_score_threshold       = 0.7; // determines how deep a prong should remove points into a track

  int nsegs = v_segs_points_preliminary.size();
  int nsegs_cm = ceil(lastPointProj-firstPointProj);
  std::vector<double> npoints_cm_seg(    nsegs_cm,0);
  std::vector<double> npoints_reject_seg(nsegs_cm,0);

  std::vector<std::vector<double>> v_prongs_COM;
  std::vector<std::vector<double>> v_prongs_dir;

  //iteratively remove points by finding regions with the most (remaining) rejected points
  while (true) {
    //Look along removed points for prongs (overlapping tracks and delta rays) and remove them
    Track track_prong;
    //re-count how many points are in each cm seg and how many points are in region around each cm seg
    int min_index = 0;
    int max_index = 0;
    for (int i=0;i<nsegs_cm;i++) { npoints_cm_seg[i]     = 0; }
    for (int i=0;i<nsegs_cm;i++) { npoints_reject_seg[i] = 0; }
    for (int i=0;i<track_reject.N;i++) { npoints_cm_seg[std::clamp( (int)floor(dot(track_reject.points[i],aAxis)-firstPointProj), 0, nsegs_cm-1 )]++; }
    for (int i=0;i<nsegs_cm;i++) {
      min_index = std::clamp((int)(i-seg_reject_length/2),0,nsegs_cm-1);
      max_index = std::clamp((int)(i+seg_reject_length/2),0,nsegs_cm-1);
      for (int j=min_index;j<max_index;j++) { npoints_reject_seg[i] += npoints_cm_seg[j]; }
    }
    //find index of most populated reject seg
    int max_reject_seg_index = std::distance(npoints_reject_seg.begin(),std::max_element(npoints_reject_seg.begin(), npoints_reject_seg.end()));
    min_index = std::clamp((int)(max_reject_seg_index-seg_reject_length/2),0,nsegs_cm-1); //relies on 1cm seg
    max_index = std::clamp((int)(max_reject_seg_index+seg_reject_length/2),0,nsegs_cm-1);
    if (npoints_reject_seg[max_reject_seg_index]<npoints_seg_reject_threshold) { break; } //if there are no segs with enough points, end the loop
    //fill potential prong with points within region along aAxis
    for (int i=0;i<track_reject.N;i++) {
      double proj = dot(track_reject.points[i],aAxis)-firstPointProj;
      if (proj>=min_index && proj<max_index) { track_prong.add_point(track_reject.points[i]); } //relies on 1cm segs, otherwise scale index by width
    }
    if(track_prong.N<npoints_seg_reject_threshold){ break;}
    Track track_prong_initial;
      for (int i=0;i<track_prong.N;i++) { track_prong_initial.add_point(track_prong.points[i]); }
    

    //use PCA to get primary direction and COM of prong
    std::vector<double> prong_COM;
    std::vector<double> prong_dir;
    std::vector<double> prong_evals; //eigenvalues to PCA fit
    prong_dir = fitPCA(track_prong,prong_COM,prong_evals)[0];
    track_prong.add_point( generate_track_point(prong_COM, averaged_COM_preliminary, v_segs_tangent, aAxis) ); //add point along ceneter of track to stabilize fit
    prong_dir = fitPCA(track_prong,prong_COM,prong_evals)[0];

    //re-fill track with any point with a perpendicular offset less than first threshold and re-fit PCA
    track_prong.clear();
    std::vector<double> v_perpU = cross(prong_dir,aAxis);
    std::vector<double> v_perpV = cross(prong_dir,v_perpU);
    v_perpU = scale(v_perpU,1./norm(v_perpU));
    v_perpV = scale(v_perpV,1./norm(v_perpV));
    std::vector<double> prong_COM_perp = { 0, dot(prong_COM,v_perpU) , dot(prong_COM,v_perpV) };
    for (int i=0;i<track_reject.N;i++) {
      std::vector<double> reject_point_perp = { 0, dot(track_reject.points[i],v_perpU) , dot(track_reject.points[i],v_perpV) };
      if (norm(diff(prong_COM_perp,reject_point_perp))<prong_width_threshold) { track_prong.add_point(track_reject.points[i]); } //relies on 1cm segs, otherwise scale index by width
    }
    

    //check whether candidate prong has enough points, otherwise ignore this candidate prong
    if (track_prong.N<npoints_seg_reject_threshold/2) {
      //always remove reject points used in PCA fit
      for (int i=0;i<track_prong_initial.N;i++) { track_reject.remove_seg(track_prong_initial.points[i],1); }
      continue;
    }
    prong_dir = fitPCA(track_prong,prong_COM,prong_evals)[0];
    track_prong.add_point( generate_track_point(prong_COM, averaged_COM_preliminary, v_segs_tangent, aAxis) ); //add point along ceneter of track to stabilize fit
    prong_dir = fitPCA(track_prong,prong_COM,prong_evals)[0];
    
    //re-fill track with points with offset within second threshold and re-fit PCA
    double prong_offset_90p = get_offset_90p(track_prong.points, prong_dir, aAxis, prong_COM, prong_COM_perp, v_perpU, v_perpV, true);
    track_prong.clear();
    for (int i=0;i<track_reject.N;i++) {
      std::vector<double> reject_point_perp = { 0, dot(track_reject.points[i],v_perpU) , dot(track_reject.points[i],v_perpV) };
      if (norm(diff(prong_COM_perp,reject_point_perp))<prong_offset_90p*prong_width_ratio_threshold) { track_prong.add_point(track_reject.points[i]); } //relies on 1cm segs, otherwise scale index by width
    }
    if (track_prong.N<npoints_seg_reject_threshold/2) {
      //always remove reject points used in PCA fit
      for (int i=0;i<track_prong_initial.N;i++) { track_reject.remove_seg(track_prong_initial.points[i],1); }
      continue;
    }
    prong_dir = fitPCA(track_prong,prong_COM,prong_evals)[0];
    track_prong.add_point( generate_track_point(prong_COM, averaged_COM_preliminary, v_segs_tangent, aAxis) ); //add point along ceneter of track to stabilize fit
    prong_dir = fitPCA(track_prong,prong_COM,prong_evals)[0];

    //find which COM point is closest to fit line (use projection onto plane)
    std::vector<double> v_seg_COM_offset;
    for (int i=0;i<nsegs;i++) {
      std::vector<double> seg_COM_perp = { 0, dot(averaged_COM_preliminary[i],v_perpU), dot(averaged_COM_preliminary[i],v_perpV) };
      v_seg_COM_offset.push_back(norm(diff(seg_COM_perp,prong_COM_perp)));
    }
    int min_seg_COM_offset_index = std::distance(v_seg_COM_offset.begin(),std::min_element(v_seg_COM_offset.begin(), v_seg_COM_offset.end()));

    //count how many prong points are on each side of COM point, if one side only has a few remove them all and re-fit PCA
    double npoints_side1 = 0;
    double npoints_side2 = 0;
    bool crossing_prong = true;
    double seg_COM_proj = dot(prong_dir, averaged_COM_preliminary[min_seg_COM_offset_index]); 
    for (int i=0;i<track_prong.N;i++) {
      double point_proj = dot(prong_dir, track_prong.points[i]);
      if ( point_proj < seg_COM_proj ) { npoints_side1++; } else { npoints_side2++; }
    }
    if (npoints_side1 < npoints_seg_reject_threshold/2) {
      crossing_prong = false;
      for (int i=track_prong.N-1;i>=0;i--) {
        double point_proj = dot(prong_dir, track_prong.points[i]);
        if ( point_proj < seg_COM_proj ) { track_prong.remove_seg_at(i,1); }
      }
    }
    
   
    if (npoints_side2 < npoints_seg_reject_threshold/2) {
      crossing_prong = false;
      for (int i=track_prong.N-1;i>=0;i--) {
        double point_proj = dot(prong_dir, track_prong.points[i]);
        if ( point_proj > seg_COM_proj ) { track_prong.remove_seg_at(i,1); }
      }
    }
    if (!crossing_prong) {
      if (track_prong.N<npoints_seg_reject_threshold/2) {
        //always remove reject points used in PCA fit
        for (int i=0;i<track_prong_initial.N;i++) { track_reject.remove_seg(track_prong_initial.points[i],1); }
        continue;
      }
      prong_dir = fitPCA(track_prong, prong_COM, prong_evals)[0];
    }

    //check whether final fit line meets thresholds for removal based on angle with tangent line and eigenvalue
    int iseg = (max_reject_seg_index+0.5)/seg_length; //relies on 1cm seg
    double prong_angle = acos(abs(dot(prong_dir,v_segs_tangent[iseg])));
    
   
    
    if ( prong_angle > prong_angle_threshold && track_prong.N >= npoints_seg_reject_threshold/2) {
    
    
      v_prongs_COM.push_back(prong_COM);
    
      v_prongs_dir.push_back(prong_dir);
      
      prong_offset_90p = 1.15*get_offset_90p(track_prong.points, prong_dir, aAxis, prong_COM, prong_COM_perp, v_perpU, v_perpV, true);

      //find which COM point is closest to fit line (use projection onto plane)
    
      v_seg_COM_offset.clear();
      for (int i=0;i<nsegs;i++) {
    
	std::vector<double> seg_COM_perp = { 0, dot(averaged_COM_preliminary[i],v_perpU), dot(averaged_COM_preliminary[i],v_perpV) };
        v_seg_COM_offset.push_back(norm(diff(seg_COM_perp,prong_COM_perp)));
      }
 
      min_seg_COM_offset_index = std::distance(v_seg_COM_offset.begin(),std::min_element(v_seg_COM_offset.begin(), v_seg_COM_offset.end()));
      double seg_COM_prong_proj   = dot(prong_dir, averaged_COM_preliminary[min_seg_COM_offset_index]); 
      double prong_COM_prong_proj = dot(prong_dir, prong_COM); 
 

      //iterate over local segs, evaluate whether each point is part of the prong or not
      min_index = std::clamp((int)(min_seg_COM_offset_index-seg_reject_length/seg_length/2),0,nsegs-1);
      max_index = std::clamp((int)(min_seg_COM_offset_index+seg_reject_length/seg_length/2),0,nsegs-1);
      for (int iseg=min_index;iseg<=max_index;iseg++) {
        std::vector<double> v_seg_perpU = cross(v_segs_tangent[iseg],aAxis);
        std::vector<double> v_seg_perpV = cross(v_segs_tangent[iseg],v_seg_perpU);
        v_seg_perpU = scale(v_seg_perpU,1./norm(v_seg_perpU));
        v_seg_perpV = scale(v_seg_perpV,1./norm(v_seg_perpV));
	std::vector<double> seg_COM_perp = { 0, dot(averaged_COM_preliminary[iseg],v_seg_perpU), dot(averaged_COM_preliminary[iseg],v_seg_perpV) };
        int npoints_seg = v_segs_points_preliminary[iseg].size();
	for (int ipoint=npoints_seg-1;ipoint>=0;ipoint--) { //iterate in reverse order because we are removing points and don't want to mess up indexing
          std::vector<double> point = v_segs_points_preliminary[iseg][ipoint];
          double point_prong_proj = dot(prong_dir, point);
          std::vector<double> v_point_prong_perp = { 0, dot(point,v_perpU),     dot(point,v_perpV)     };
          std::vector<double> v_point_seg_perp   = { 0, dot(point,v_seg_perpU), dot(point,v_seg_perpV) };
          double offset_prong = norm(diff(v_point_prong_perp,prong_COM_perp));
	  double offset_seg   = norm(diff(v_point_seg_perp,  seg_COM_perp));
          double prong_score = -1*std::pow(offset_prong / prong_offset_90p,       4);
	  double seg_score   = -1*std::pow(offset_seg   / v_segs_offset_90p[iseg],4);
          bool point_prong_sameside = ((prong_COM_prong_proj-seg_COM_prong_proj) * (point_prong_proj-seg_COM_prong_proj)) > 0;
	  if (prong_evals[0] > 0.7 && prong_score > seg_score + prong_score_threshold && (crossing_prong || point_prong_sameside)) {
            v_segs_points_preliminary[iseg].erase(v_segs_points_preliminary[iseg].begin()+ipoint);
	  }
	} //loop points
      } //loop segs
   
    } //if removing main track points

    //always remove reject points used in PCA fit (don't remove last point that was artificially added along center of track)
    for (int i=0;i<track_prong.N-1;i++) { track_reject.remove_seg(track_prong.points[i],1); }
     
  } //while loop removing prongs
 // std::cout<<"v_prongs_COM.size ,v_prongs_dir.size "<<v_prongs_COM.size()<<", "<<v_prongs_dir.size()<<std::endl;
  
  return { v_prongs_COM, v_prongs_dir };
}


//Takes in a vector of points as well as muon start,end reco coordinates and returns a refined subset of those points that belong to the muon
//aimed at removing delta rays and crossing tracks
std::vector<std::vector<std::vector<double>>> select_points_muon_refined (std::vector<std::vector<double>> v_points_all, std::vector<std::vector<double>> v_points_candidate, std::vector<double> muon_start, std::vector<double> muon_end) {

  double seg_length           = 1.0;
  double local_range_v1       = 7;
  double local_range_v2       = 7;
  double local_range_v3       = 5;
  double local_range_v4       = 5;
  double trim_radius_v1       = 1.8; //1.5
  double trim_radius_v2       = 1.6; //1.3
  double trim_offset_ratio_v3 = 1.25;
  double trim_offset_ratio_v4 = 1.1;
  int local_range_index_v1 = (int)(local_range_v1 / seg_length);
  int local_range_index_v2 = (int)(local_range_v2 / seg_length);
  int local_range_index_v3 = (int)(local_range_v3 / seg_length);
  int local_range_index_v4 = (int)(local_range_v4 / seg_length);
	

  //fill track
  Track track;
  int npoints_all = v_points_all.size();
  for(int i=0; i<npoints_all; i++){ track.add_point(v_points_all[i]); }

  //get PCA axes for track
  std::vector<double> com         = {0,0,0};
  std::vector<double> eigenvalues = {0,0,0};
  std::vector<std::vector<double>> axes = fitPCA(track, com, eigenvalues);
  std::vector<double> aAxis = axes[0];
  std::vector<double> bAxis = axes[1];
  std::vector<double> cAxis = axes[2];
  //flip axis to be along muon direction
  std::vector<double> vec_muon = { muon_end[0]-muon_start[0], muon_end[1]-muon_start[1], muon_end[2]-muon_start[2] };
  if(dot(vec_muon,aAxis) < 0){ aAxis = scale(aAxis,-1.); }
  sort_points(track.points, aAxis);
  double firstPointProj = dot(track.points.front(),aAxis);
  double lastPointProj  = dot(track.points.back(),aAxis);
  std::vector<double> first_point = track.points.front();
  track.clear();

  int nsegs = ceil((lastPointProj-firstPointProj)/seg_length);
  std::vector<std::vector<std::vector<double>>> v_segs_points; //vector of segs of xyz points
  std::vector<std::vector<double>> v_segs_tangent;    //vector of tangent lines for each segment
  std::vector<double> v_segs_offset_90p; //vector of 90% offset vals for each segment
  for (int i=0;i<nsegs;i++) {
    v_segs_points.push_back({});
    v_segs_tangent.push_back({});
    v_segs_offset_90p.push_back({});
  }

  //subdivide into segments of specified length
  //use rigid length, allowing empty segs over gaps
  for (int i=0;i<npoints_all;i++) {
    int iseg = (int)((dot(v_points_all[i],aAxis)-firstPointProj)/seg_length);
    v_segs_points[iseg].push_back(v_points_all[i]);
  }

  //compute COM points along track by averaging within +-7cm then trim point cloud points outside threshold of any COM points
  std::vector<std::vector<double>> averaged_COM_v1                        = compute_averaged_COM(  v_segs_points,            v_segs_tangent, aAxis, first_point, seg_length, local_range_index_v1);
  std::vector<std::vector<std::vector<double>>> v_segs_points_trimmed_v1  = trim_seg_points_radius(v_segs_points, averaged_COM_v1, trim_radius_v1);

  //repeat with second pass and smaller threshold
  std::vector<std::vector<double>> averaged_COM_v2                        = compute_averaged_COM(  v_segs_points_trimmed_v1, v_segs_tangent, aAxis, first_point, seg_length, local_range_index_v2);
  std::vector<std::vector<std::vector<double>>> v_segs_points_trimmed_v2  = trim_seg_points_radius(v_segs_points_trimmed_v1, averaged_COM_v2, trim_radius_v2);

  //repeat with third pass: smaller seg range and trim point cloud based on offset ratio from tangent lines
  std::vector<std::vector<double>> averaged_COM_v3                        = compute_averaged_COM(  v_segs_points_trimmed_v2, v_segs_tangent, aAxis, first_point, seg_length, local_range_index_v3);
  std::vector<std::vector<std::vector<double>>> v_segs_points_trimmed_v3  = trim_seg_points_offset(v_segs_points_trimmed_v2, averaged_COM_v3, v_segs_tangent, v_segs_offset_90p, local_range_index_v3, trim_offset_ratio_v3);

  //repeat with fourth pass (similar to 3rd) and smaller threshold
  std::vector<std::vector<double>> averaged_COM_v4                        = compute_averaged_COM(  v_segs_points_trimmed_v3, v_segs_tangent, aAxis, first_point, seg_length, local_range_index_v4);
  std::vector<std::vector<std::vector<double>>> v_segs_points_trimmed_v4  = trim_seg_points_offset(v_segs_points_trimmed_v3, averaged_COM_v4, v_segs_tangent, v_segs_offset_90p, local_range_index_v4, trim_offset_ratio_v4);

  //get COM and tangent lines for preliminary selection.  Call trim_seg_points_offset just to update v_segs_offset_90p and generate v_segs_points_preliminary
  std::vector<std::vector<double>> averaged_COM_preliminary               = compute_averaged_COM(  v_segs_points_trimmed_v4, v_segs_tangent, aAxis, first_point, seg_length, local_range_index_v4);
  std::vector<std::vector<std::vector<double>>> v_segs_points_preliminary = trim_seg_points_offset(v_segs_points_trimmed_v4, averaged_COM_v4, v_segs_tangent, v_segs_offset_90p, local_range_index_v4, 1e6);

  v_segs_points_trimmed_v1.clear();
  v_segs_points_trimmed_v2.clear();
  v_segs_points_trimmed_v3.clear();
  v_segs_points_trimmed_v4.clear();

  //combine segs into preliminary point cloud
  std::vector<std::vector<double>> v_points_preliminary;
  for (int iseg=0;iseg<nsegs;iseg++) {
    int npoints = v_segs_points_preliminary[iseg].size();
    for (int j=0;j<npoints;j++) {
      v_points_preliminary.push_back(v_segs_points_preliminary[iseg][j]);
      track.add_point(v_segs_points_preliminary[iseg][j]);
    }
  }
  int npoints_preliminary = v_points_preliminary.size();


  //re-fit PCA
  axes = fitPCA(track, com, eigenvalues);
  aAxis = axes[0];
  bAxis = axes[1];
  cAxis = axes[2];
  if(dot(vec_muon,aAxis) < 0){ aAxis = scale(aAxis,-1.); }
  sort_points(track.points, aAxis);
  firstPointProj = dot(track.points.front(),aAxis);
  lastPointProj  = dot(track.points.back(),aAxis);
  track.clear();


  //fill reject points track
  Track track_reject;
  int npoints_candidate = v_points_candidate.size();
  for(int i=0; i<npoints_candidate; i++) { track_reject.add_point(v_points_candidate[i]); }
  for(int i=0;i<npoints_preliminary;i++) { track_reject.remove_seg(v_points_preliminary[i],1); }

  std::vector<std::vector<std::vector<double>>> v_prongs = remove_prongs(track_reject, v_segs_points_preliminary, averaged_COM_preliminary, v_segs_tangent, v_segs_offset_90p, aAxis, seg_length, firstPointProj, lastPointProj, bAxis, cAxis);

  std::vector<std::vector<double>> v_prongs_COM = v_prongs[0];
  std::vector<std::vector<double>> v_prongs_dir = v_prongs[1];

  //form final point cloud of selected points
  std::vector<std::vector<double>> v_points_final;
  for (int i=0;i<nsegs;i++) {
    int npoints = v_segs_points_preliminary[i].size();
    for (int j=0;j<npoints;j++) {
      v_points_final.push_back(v_segs_points_preliminary[i][j]);
    }
  }

  return { v_points_preliminary, v_points_final, v_prongs_COM, v_prongs_dir };
}
*/



