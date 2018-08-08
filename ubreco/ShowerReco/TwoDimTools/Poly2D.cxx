#ifndef TWODIM_POLY2D_CXX
#define TWODIM_POLY2D_CXX

#include "Poly2D.h"

namespace twodimtools {

  Poly2D::Poly2D()
  { 
    vertices.clear(); 
    // get detector specific properties
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    _wire2cm = geom->WirePitch(0,1,0);
    _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );
    _trigoff = detp->TriggerOffset();
  }

  Poly2D::Poly2D(const std::vector< art::Ptr<recob::Hit> >& hit_v)
    : Poly2D()
  {
    SelectPolygonHitList(hit_v,vertices,0.95);
  }
  
  //------------------------------------------------
  // returns slope of line uniting points p1 and p2
  float FindSlope( const std::pair<float, float> &p1,
		   const std::pair<float, float> &p2 )
  {
    float slope = (p2.second - p1.second) / (p2.first - p1.first);
    return slope;
  }
  
  //-------------------------------------------------------------------------
  // determines if 3 points are arranged in a clock-wire oder or not
  bool Clockwise(double Ax, double Ay, double Bx, double By, double Cx, double Cy)
  {
    return (Cy - Ay) * (Bx - Ax) > (By - Ay) * (Cx - Ax);
  }
  
  //------------------------------------------------------------
  // determine if two segments intersect
  bool SegmentOverlap(double Ax, double Ay, double Bx, double By,
		      double Cx, double Cy, double Dx, double Dy)
  {
    
    bool overlap = ( (Clockwise(Ax, Ay, Cx, Cy, Dx, Dy) != Clockwise(Bx, By, Cx, Cy, Dx, Dy))
		     and (Clockwise(Ax, Ay, Bx, By, Cx, Cy) != Clockwise(Ax, Ay, Bx, By, Dx, Dy)) );
    return overlap;
  }
  
  //---------------------------------------------------------------------------------
  // return intersection point for two segments
  std::pair<float, float> GetIntersection(double Ax, double Ay, double Bx, double By,
					  double Cx, double Cy, double Dx, double Dy)
  {
    
    //get equations for two lines
    // [Ax,Ay]<--->[Bx,By] : y = s1*x+c1
    // [Cx,Cy]<--->[Dx,Dy] : y = s2*x+c2
    double s1 = (By - Ay) / (Bx - Ax);
    double s2 = (Dy - Cy) / (Dx - Cx);
    double c1 = By - s1 * Bx;
    double c2 = Dy - s2 * Dx;
    
    double Xintersection = (c2 - c1) / (s2 - s1);
    double Yintersection = s1 * Xintersection + c1;
    std::pair<float, float> intersection;
    intersection = std::make_pair(Xintersection, Yintersection);
    
    return intersection;
  }
  
  //------------------------------------------------------------------
  // Poly2D constructor returning intersection of 2 polygons
  Poly2D::Poly2D(const Poly2D &poly1, const Poly2D &poly2)
  {
    
    //figure out if the two polygons overlap at all
    if ( !(poly1.PolyOverlap(poly2)) ) {
      std::vector< std::pair<float, float> > nullpoint;
      vertices = nullpoint;
      return;
    }
    
    //The overlap polygon is made up by:
    //1) all points of poly1 in poly2
    //2) all points of poly2 in poly1
    //3) all intersection points between segments
    
    //make a new set of points and add points
    //as listed above, if found.
    std::vector<std::pair<float, float> > IntersectionPoints;
    //1)
    for (unsigned int p1 = 0; p1 < poly1.Size(); p1++) {
      if ( poly2.PointInside( poly1.Point(p1) ) ) { IntersectionPoints.push_back( poly1.Point(p1) ); }
    }
    //2)
    for (unsigned int p2 = 0; p2 < poly2.Size(); p2++) {
      if ( poly1.PointInside( poly2.Point(p2) ) ) { IntersectionPoints.push_back( poly2.Point(p2) ); }
    }
    //3)
    //FIND SEGMENT INTERSECTIONS
    for (unsigned int i = 0; i < poly1.Size(); i++) {
      for (unsigned int j = 0; j < poly2.Size(); j++) {
	if (SegmentOverlap( poly1.Point(i).first, poly1.Point(i).second,
			    poly1.Point(i + 1).first, poly1.Point(i + 1).second,
			    poly2.Point(j).first, poly2.Point(j).second,
			    poly2.Point(j + 1).first, poly2.Point(j + 1).second) )
	  //segments overlap...add intersection to list
	  IntersectionPoints.push_back( GetIntersection( poly1.Point(i).first, poly1.Point(i).second,
							 poly1.Point(i + 1).first, poly1.Point(i + 1).second,
							 poly2.Point(j).first, poly2.Point(j).second,
							 poly2.Point(j + 1).first, poly2.Point(j + 1).second) );
      }//for all segments in poly2
    }//for all segments in poly1
    
    vertices = IntersectionPoints;
    return;
  }
  
  //---------------------------
  float Poly2D::Area() const
  {
    //how? here:
    //http://www.mathsisfun.com/geometry/area-irregular-polygons.html
    float area = 0;
    for (unsigned int i = 0; i < vertices.size(); i++) {
      if ( i < (vertices.size() - 1) )
	area += (((vertices.at(i)).second) + ((vertices.at(i + 1)).second)) * (((vertices.at(i)).first) - ((vertices.at(i + 1)).first)) * 0.5;
      if ( i == (vertices.size() - 1) )
	area += (((vertices.at(i)).second) + ((vertices.at(0)).second)) * (((vertices.at(i)).first) - ((vertices.at(0)).first)) * 0.5;
    }
    if (area < 0.0)
      area = -area;
    return area;
  }
  
  //--------------------------------
  float Poly2D::Perimeter() const
  {
    
    float perimeter = 0.;
    
    for (unsigned int i = 0; i < vertices.size(); i++) {
      if ( i < (vertices.size() - 1) )
	perimeter += ( (vertices.at(i).second - vertices.at(i + 1).second) *
		       (vertices.at(i).second - vertices.at(i + 1).second) +
		       (vertices.at(i).first - vertices.at(i + 1).first) *
		       (vertices.at(i).first - vertices.at(i + 1).first) );
      if ( i == (vertices.size() - 1) )
	perimeter += ( (vertices.at(i).second - vertices.at(0).second) *
		       (vertices.at(i).second - vertices.at(0).second) +
		       (vertices.at(i).first - vertices.at(0).first) *
		       (vertices.at(i).first - vertices.at(0).first) );
    }
    
    return sqrt(perimeter);
  }
  
  //------------------------------------------------------------------
  const std::pair<float, float>& Poly2D::Point(unsigned int p) const
  {
    //This function returns the vertex under consideration
    //as a std::pair<float,float> Returns vertex for argument
    //from 0 to N-1. For input N = number of sides then
    //the first vertex is returned
    if (p < vertices.size())
      return vertices.at(p);
    else if (p == vertices.size())
      return vertices.at(0);
    else {
      std::cout << "Out of bounds of Polygon!" << std::endl;
      return vertices.at(0);
    }
    
  }
  
  //------------------------------------------------------------------------
  // apply translation and rotation to a polygon
  std::pair<float, float> Poly2D::Project(const std::pair<float, float> &p,
					     float theta) const
  {
    
    std::pair<float, float> range(10000, 0);
    std::pair<float, float> ptmp;
    
    for (unsigned int i = 0; i < vertices.size(); i++) {
      //Translation
      //translating each vertex so that origin is on first vertex on polygon's edge being considered
      ptmp = std::make_pair(   (vertices.at(i)).first - p.first  ,   (vertices.at(i)).second - p.second   );
      //Rotation
      //instead of rotating each (x,y) edge (slow) just find nex x-position which gives us information
      //on the projection of that vertex on the line we are considering
      // +x direction is from vertex in consideration (vertex 'i' in loop) to next vertex
      //now find the x-coordinate of that vertex after it is rotated such that edge is now + x axis
      float xnew = (ptmp.first) * cos(theta) + (ptmp.second) * sin(theta);
      //finally calculate range of projection on x-axis: look at every x position and compare it to range
      if ( xnew < range.first )
	range.first = xnew;
      if ( xnew > range.second )
	range.second = xnew;
    }
    return range;
  }

  //---------------------------------------------------------------
  bool Poly2D::Overlap(const Poly2D& poly2) const
    
  {
    
    bool overlap = false;

    
    for (size_t j=0; j < poly2.Size(); j++){
      
      if ( this->PointInside(poly2.Point(j)) == true) {
	overlap = true;
	break;
      }
    }

    return overlap;
  }
  
  //---------------------------------------------------------------
  bool Poly2D::Overlap(float slope,
		       const Poly2D &poly2,
		       const std::pair<float, float> &origin) const
  {
    //translate and rotate both polygons
    float theta = tan(slope);
    //here we translate origin, rotate and find x-coordinates and find range of projection on edge line
    std::pair<float, float> range1 = this->Project(origin, theta);
    std::pair<float, float> range2 = poly2.Project(origin, theta);
    //std::cout << "range 1: " << range1.first << " " << range1.second << std::endl;
    //std::cout << "range 2: " << range2.first << " " << range2.second << std::endl;
    //if the two projections don't overlap --> no overlap!
    if ( ( ((range1.first <= range2.second) and ( range1.first >= range2.first )) or ((range1.second <= range2.second) and ( range1.second >= range2.first )) ) or ( ((range2.first <= range1.second) and ( range2.first >= range1.first )) or ((range2.second <= range1.second) and ( range2.second >= range1.first )) ) )
      return true;     //yes...they overlap
    else
      return false;    //no....they do not overlap
  }
  
  //-------------------------------------------------------
  bool Poly2D::PolyOverlap(const Poly2D &poly2) const
  {
    
    //start from first pair in vector then check all edges.
    //edges are (0,1), (1,2), ..., (N,N-1) each pair a pair
    //of vertexes
    for (unsigned int i = 0; i < this->Size(); i++) { //loop over first polygon's vertices
      //find projection line's slope
      //line: y=ax+b --- slope is a variable
      float slope;
      slope = FindSlope( this->Point(i) , this->Point(i + 1) );
      //if there is even one no-overlap
      //need to exit and return no overlap!
      if (! (this->Overlap( slope, poly2, this->Point(i) )) )
	return false;
    }//loop over first polygon vertices
    
    //do the exact same thing but reversing polygon role
    for (unsigned int i = 0; i < poly2.Size(); i++) { //loop over first polygon
      float slope;
      slope = FindSlope( poly2.Point(i) , poly2.Point(i + 1) );
      if (!(poly2.Overlap( slope, *this, poly2.Point(i) )) )
	return false;
    }//loop over second polygon vertices
    return true;
  }
  
  //---------------------------------------------------------------
  bool Poly2D::PolyOverlapSegments(const Poly2D &poly2) const
  {
    //if contained in one another then they also overlap:
    if ( (this->Contained(poly2)) or (poly2.Contained(*this)) ) {
      return true;
    }
    //loop over the two polygons checking wehther
    //two segments ever intersect
    for (unsigned int i = 0; i < this->Size(); i++) {
      for (unsigned int j = 0; j < poly2.Size(); j++) {
	if (SegmentOverlap( this->Point(i).first, this->Point(i).second,
			    this->Point(i + 1).first, this->Point(i + 1).second,
			    poly2.Point(j).first, poly2.Point(j).second,
			    poly2.Point(j + 1).first, poly2.Point(j + 1).second) ) {
	  return true;
	}
      }
    }
    return false;
  }
  
  //--------------------------------------------------------------------
  bool Poly2D::PointInside(const std::pair<float, float> &point) const
  {
    
    //any ray originating at point will cross polygon
    //even number of times if point outside
    //odd number of times if point inside
    int intersections = 0;
    for (unsigned int i = 0; i < this->Size(); i++) {
      if ( SegmentOverlap( this->Point(i).first, this->Point(i).second,
			   this->Point(i + 1).first, this->Point(i + 1).second,
			   10000.0, 10000.0,
			   point.first, point.second) )
	intersections += 1;
    }
    if ( (intersections % 2) == 0 )
      return false;
    else
      return true;
    
  }
  
  //-----------------------------------------------------
  bool Poly2D::Contained(const Poly2D &poly2) const
  {
    
    //loop over poly2 checking wehther
    //points of poly2 all inside poly1
    for (unsigned int i = 0; i < poly2.Size(); i++) {
      if ( !(this->PointInside( poly2.Point(i)) ) )
	return false;
    }
    
    return true;
    
  }
  
  //-------------------------------
  void Poly2D::UntanglePolygon()
  {
    
    //loop over edges
    for ( unsigned int i = 0; i < vertices.size() - 1; i++) {
      double Ax = vertices.at(i).first;
      double Ay = vertices.at(i).second;
      double Bx = vertices.at(i + 1).first;
      double By = vertices.at(i + 1).second;
      //loop over edges that have not been checked yet
      for (unsigned int j = i + 2; j < vertices.size() - 1; j++) {
	//avoid consecutive segments
	if ( vertices.at(i) == vertices.at(j + 1) )
	  continue;
	else {
	  double Cx = vertices.at(j).first;
	  double Cy = vertices.at(j).second;
	  double Dx = vertices.at(j + 1).first;
	  double Dy = vertices.at(j + 1).second;
	  
	  if ( SegmentOverlap( Ax, Ay, Bx, By, Cx, Cy, Dx, Dy ) ) {
	    std::pair<float, float> tmp = vertices.at(i + 1);
	    vertices.at(i + 1) = vertices.at(j);
	    vertices.at(j) = tmp;
	    //swapped polygon, now do recursion to make sure
	    this->UntanglePolygon();
	  }//if crossing
	}
      }//second loop
    }//first loop
    
  }

  /*  
  ///Calculate the opening angle at the specified vertex:
  float Poly2D::InteriorAngle(unsigned int p) const {
    // Get the specified point and the two adjacent points
    if (p > vertices.size()) {
      return -9999999;
    }
    
    // std::cout << "Getting angle of point " << p << " at (" << vertices.at(p).first
    //           << ", " << vertices.at(p).second << ")" << std::endl;
    
    auto geoHelper = larutil::GeometryHelper::GetME();
    
    // Need at least 3 points:
    if (vertices.size() < 3)
      return 0.0;
    else {
      unsigned int next, prev;
      // Check if this is the very last point:
      next = (p + 1) % Size();
      prev = (p - 1 + Size()) % Size();
      // Now actually calculate the cosine of the angle with the dot product:
      larutil::Point2D point_prev, point_next, point_p;
      point_prev.w       = Point(prev).first; point_prev.t = Point(prev).second;
      point_next.w       = Point(next).first; point_next.t = Point(next).second;
      point_p.w          = Point(p).first;    point_p.t = Point(p).second;
      
      
      // Determine if this is an interior or exterior angle:
      if (vertices.size() == 3) {
	// All triangle angles are interior:
	return acos(geoHelper -> GetCosAngleBetweenLines(point_p, point_next, point_prev));
      }
      else {
	// Get the next point in the polygon
	std::vector<int> otherPoints;
	otherPoints.reserve(Size());
	for (int point = 0; point < Size(); point ++) {
	  if (point != p) {
	    otherPoints.push_back(point);
	  }
	}
	
	
	
	// std::cout << "\tPrev:   (" << vertices.at(prev).first << ", " << vertices.at(prev).second << ")" << std::endl;
	// std::cout << "\tNext:   (" << vertices.at(next).first << ", " << vertices.at(next).second << ")" << std::endl;
	
	// Determine if convex or concave by figuring out if it is inside or outside the sub-polygon without it
	// Basically using the algorithm from here http://stackoverflow.com/questions/217578/point-in-polygon-aka-hit-test
	// Need to know this in order to correctly calculate interior angle
	
	
	
	int hit_count = 0;
	for (unsigned int i = 0; i < otherPoints.size(); i++) {
	  int point = otherPoints[i];
	  int point2;
	  if (point == otherPoints.back()) {
	    point2 = otherPoints.front();
	  }
	  else {
	    point2 = otherPoints[i + 1];
	  }
	  if (SegmentOverlap(vertices.at(p).first, vertices.at(p).second,
			     -10, -10, //Compare it to an unphysical point to ensure it's not in the polygon
			     vertices.at(point).first,
			     vertices.at(point).second,
			     vertices.at(point2).first,
			     vertices.at(point2).second))
	    {
	      hit_count ++;
	    }
	}
	
	if (hit_count % 2 == 0) {
	  // Then it crossed an even number of points and is therefore not in the subpolygon.  Convex
	  // std::cout << "\tConvex, cosAngle is "
	  //           << geoHelper -> GetCosAngleBetweenLines(point_p, point_next, point_prev)
	  //           << ", returning "
	  //           <<  acos(geoHelper -> GetCosAngleBetweenLines(point_p, point_next, point_prev)) << std::endl;
	  return acos(geoHelper -> GetCosAngleBetweenLines(point_p, point_next, point_prev));
	}
	else
	  {
	    // std::cout << "\tConcave, cosAngle is  "
	    //           << geoHelper -> GetCosAngleBetweenLines(point_p, point_next, point_prev)
	    //           << ", returning "
	    //           <<  2*M_PI - acos(geoHelper -> GetCosAngleBetweenLines(point_p, point_next, point_prev)) << std::endl;
	    return 2 * M_PI - acos(geoHelper -> GetCosAngleBetweenLines(point_p, point_next, point_prev));
	  }
      }
    }
  }
  
  bool operator==(const Poly2D& lhs, const Poly2D& rhs) {
    return lhs.vertices == rhs.vertices;
  }
  
  bool operator!=(const Poly2D& lhs, const Poly2D& rhs) {
    return !(lhs == rhs);
  }
  
  std::ostream &operator<<(std::ostream &out, Poly2D poly) {
    
    if (poly.vertices.size( ) == 0) {
      out << "(NULL POLYGON)";
      return out;
    }
    
    for (unsigned int i = 0; i < poly.vertices.size(); ++i)
      {
	out << "(" << poly.vertices.at(i).first << ", " << poly.vertices.at(i).second << ")";
	if (i != poly.vertices.size() - 1 )
	  out << ", ";
      }
    
    return out;
  }
  */
  
  void Poly2D::SelectPolygonHitList(const std::vector< art::Ptr<recob::Hit> >& hit_v,
				    std::vector< std::pair<float,float> > &edgehits,
				    double frac) const
  {
    
    // if hit list is empty get out of here!
    if (!(hit_v.size())) {
      std::cout << "Provided empty hit list!" << std::endl;
      return;
    }
    
    // if the fraction is > 1 then use 1...should not be larger
    // frac is the fraction of charge in the hit list
    // than needs to be included in the Polygon
    if (frac > 1) { frac = 1; }
    
    // clear list of hits that define the edges of the polygon
    edgehits.clear();
    
    // Define subset of hits to define polygon
    std::map<double, size_t > hitmap;
    
    // define a parameter that stores the total charge in the cluster
    double qtotal = 0;
    for (size_t i=0; i < hit_v.size(); i++) {
      hitmap.insert(std::pair<double, size_t>(hit_v.at(i)->Integral(), i));
      qtotal += hit_v.at(i)->Integral();
    }
    // define a parameter to store the charge that will be within the polygon
    double qintegral = 0;
    std::vector<std::pair<float,float>> ordered_hits;
    ordered_hits.reserve(hit_v.size());
    for (auto hiter = hitmap.rbegin(); qintegral <= qtotal * frac && hiter != hitmap.rend(); ++hiter) {
      qintegral += (*hiter).first;
      auto hit = hit_v.at( (*hiter).second );
      float w = (float)(hit->WireID().Wire * _wire2cm);
      float t = (float)((hit->PeakTime() - _trigoff) * _time2cm + hit->WireID().Plane * 0.3);
      auto hitpt = std::make_pair(w, t);
      ordered_hits.push_back( hitpt );
    }
    
    // Define container to hold found polygon corner PxHit index & distance
    std::vector<size_t> hit_index(8, 0);
    std::vector<double> hit_distance(8, 1e9);
    
    // Loop over hits and find corner points in the plane view
    // Also fill corner edge points
    std::vector< std::pair<float,float> > edges(4, std::make_pair(0,0) );
    double wire_max = 10000.;
    double time_max = 10000.;
    
    for (size_t index = 0; index < ordered_hits.size(); ++index) {
      
      double dist = 0;
      
      auto hitw = ordered_hits.at(index).first;
      auto hitt = ordered_hits.at(index).second;
      
      // First thing to do:
      // Find the hits that have the largest/smallest wire number and time
      // these will define the first (up to) 4 boundaries of our polygon
      
      // Comparison w/ (Wire,0)
      dist = hitt;
      if (dist < hit_distance.at(1)) {
	hit_distance.at(1) = dist;
	hit_index.at(1) = index;
	edges.at(0).second = hitt;
	edges.at(1).second = hitt;
      }
      
      // Comparison w/ (WireMax,Time)
      dist = wire_max - hitw;
      if (dist < hit_distance.at(3)) {
	hit_distance.at(3) = dist;
	hit_index.at(3) = index;
	edges.at(1).first = hitw;
	edges.at(2).first = hitw;
      }
      
      // Comparison w/ (Wire,TimeMax)
      dist = time_max - hitt;
      if (dist < hit_distance.at(5)) {
	hit_distance.at(5) = dist;
	hit_index.at(5) = index;
	edges.at(2).second = hitt;
	edges.at(3).second = hitt;
      }
      
      // Comparison w/ (0,Time)
      dist = hitw;
      if (dist < hit_distance.at(7)) {
	hit_distance.at(7) = dist;
	hit_index.at(7) = index;
	edges.at(0).first = hitw;
	edges.at(3).first = hitw;
      }
    }
    
    // next find the hits that are closest to the 3 corners of the rectangle
    for (size_t index = 0; index < ordered_hits.size(); ++index) {
      
      auto hitw = ordered_hits.at(index).first;
      auto hitt = ordered_hits.at(index).second;
      
      double dist = 0;
      // Comparison w/ (0,0)
      dist = pow((hitt - edges.at(0).second), 2) + pow((hitw - edges.at(0).first), 2);
      if (dist < hit_distance.at(0)) {
	hit_distance.at(0) = dist;
	hit_index.at(0) = index;
      }
      
      // Comparison w/ (WireMax,0)
      dist = pow((hitt - edges.at(1).second), 2) + pow((hitw - edges.at(1).first), 2);
      if (dist < hit_distance.at(2)) {
	hit_distance.at(2) = dist;
	hit_index.at(2) = index;
      }
      
      // Comparison w/ (WireMax,TimeMax)
      dist = pow((hitt - edges.at(2).second), 2) + pow((hitw - edges.at(2).first), 2);
      if (dist < hit_distance.at(4)) {
	hit_distance.at(4) = dist;
	hit_index.at(4) = index;
      }
      
      // Comparison w/ (0,TimeMax)
      dist = pow((hitt - edges.at(3).second), 2) + pow((hitw - edges.at(3).first), 2);
      if (dist < hit_distance.at(6)) {
	hit_distance.at(6) = dist;
	hit_index.at(6) = index;
      }
      
    }
    // Loop over the resulting hit indexes and append unique hits to define the polygon to the return hit list
    std::set<size_t> unique_index;
    std::vector<size_t> candidate_polygon;
    candidate_polygon.reserve(9);
    //    std::cout << "Original polygon: " << std::endl;
    for (auto &index : hit_index) {
      
      if (unique_index.find(index) == unique_index.end()) {
	unique_index.insert(index);
	candidate_polygon.push_back(index);
      }
    }
    for (auto &index : hit_index) {
      candidate_polygon.push_back(index);
      break;
    }
    
    // we should only have a maximum of 8 edges for the polygon!
    if (unique_index.size() > 8) {
      std::cout << "Size of the polygon > 8!" << std::endl;
      return;
    }
    
    //Untangle Polygon
    for (auto const& polyhit : candidate_polygon) {
      if (polyhit >= ordered_hits.size()) {
	std::cout << "out of range! quit..." << std::endl;
	return;
      }
    }
    if ( (ordered_hits.size() < 2) || (candidate_polygon.size() < 2) )
      return;
    candidate_polygon = OrderPolygonEdges( ordered_hits, candidate_polygon);

    
    edgehits.clear();
    for ( unsigned int i = 0; i < (candidate_polygon.size() - 1); i++) {
      edgehits.push_back(ordered_hits.at(candidate_polygon.at(i)));
    }
    
    //check that polygon does not have more than 8 sides
    if (unique_index.size() > 8) {
      std::cout << "Size of the polygon > 8!" << std::endl;
      return;
    }
    
    return;
  }
  
  
  std::vector<size_t>  Poly2D::OrderPolygonEdges(std::vector<std::pair<float,float>> ordered_hits,
						 std::vector<size_t> candidate_polygon) const
{
  
  //loop over edges
  for ( unsigned int i = 0; i < (candidate_polygon.size() - 1); i++) {
    double Ax = ordered_hits.at(candidate_polygon.at(i)).first;
    double Ay = ordered_hits.at(candidate_polygon.at(i)).second;
    double Bx = ordered_hits.at(candidate_polygon.at(i + 1)).first;
    double By = ordered_hits.at(candidate_polygon.at(i + 1)).second;
    //loop over edges that have not been checked yet...
    //only ones furhter down in polygon
    for ( unsigned int j = i + 2; j < (candidate_polygon.size() - 1); j++) {
      //avoid consecutive segments:
      if ( candidate_polygon.at(i) == candidate_polygon.at(j + 1) )
        continue;
      else {
        double Cx = ordered_hits.at(candidate_polygon.at(j)).first;
        double Cy = ordered_hits.at(candidate_polygon.at(j)).second;
        double Dx = ordered_hits.at(candidate_polygon.at(j + 1)).first;
        double Dy = ordered_hits.at(candidate_polygon.at(j + 1)).second;

        if ( (Clockwise(Ax, Ay, Cx, Cy, Dx, Dy) != Clockwise(Bx, By, Cx, Cy, Dx, Dy))
             and (Clockwise(Ax, Ay, Bx, By, Cx, Cy) != Clockwise(Ax, Ay, Bx, By, Dx, Dy)) ) {
          size_t tmp = candidate_polygon.at(i + 1);
          candidate_polygon.at(i + 1) = candidate_polygon.at(j);
          candidate_polygon.at(j) = tmp;
          //check that last element is still first (to close circle...)
          candidate_polygon.at(candidate_polygon.size() - 1) = candidate_polygon.at(0);
          //swapped polygon...now do recursion to make sure
          return OrderPolygonEdges( ordered_hits, candidate_polygon);
        }//if crossing
      }
    }//second loop
  }//first loop
  return candidate_polygon;
}

bool Poly2D::Clockwise(const double& Ax, const double& Ay,
		       const double& Bx, const double& By,
		       const double& Cx, const double& Cy) const
{

  return (Cy - Ay) * (Bx - Ax) > (By - Ay) * (Cx - Ax);
}


}// namespace

#endif
