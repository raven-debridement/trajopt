#include "geometry_2d.hpp"

namespace Geometry2D {
  
  inline bool feq(double a, double b) {
    return fabs(a - b) < eps;
  }

  inline bool fleq(double a, double b) {
    return a - b < eps;
  }

  inline bool fzero(double a) {
    return feq(a, 0);
  }

  inline bool fnpos(double a) {
    return fleq(a, 0);
  }

  inline bool fnneg(double a) {
    return fleq(0, a);
  }

  inline double direction(const Point& pi, const Point& pj, const Point& pk) {
    return (pk.x() - pi.x()) * (pj.y() - pi.y()) -  (pk.y() - pi.y()) * (pj.x() - pi.x());
  }

  inline bool on_segment(const Point& pi, const Point& pj, const Point& pk) {
    return (fleq(fmin(pi.x(), pj.x()), pk.x()) && fleq(pk.x(), fmax(pi.x(), pj.x())) &&
            fleq(fmin(pi.y(), pj.y()), pk.y()) && fleq(pk.y(), fmax(pi.y(), pj.y())));
  }

  inline double dist(const Point& a, const Point& b) {
    return (a - b).norm();
  }

  double dist_to_line(const Point& p, const Point& v, const Point& w) {
    // Return minimum distance between line segment vw and point p
    double l2 = (v - w).squaredNorm();  // i.e. |w-v|^2 -  avoid a sqrt
    if (l2 == 0.0) return dist(p, v);
    // Consider the line extending the segment, parameterized as v + t (w - v).
    // We find projection of point p onto the line. 
    // It falls where t = [(p-v) . (w-v)] / |w-v|^2
    double t = (p - v).dot(w - v) / l2;
    Point projection = v + t * (w - v);
    return dist(p, projection);
  }

  double dist_to_segment(const Point& p, const Point& v, const Point& w) {
    // Return minimum distance between line segment vw and point p
    double l2 = (v - w).squaredNorm();  // i.e. |w-v|^2 -  avoid a sqrt
    if (l2 == 0.0) return dist(p, v);
    // Consider the line extending the segment, parameterized as v + t (w - v).
    // We find projection of point p onto the line. 
    // It falls where t = [(p-v) . (w-v)] / |w-v|^2
    double t = (p - v).dot(w - v) / l2;
    if (t < 0.0) return dist(p, v);
    else if (t > 1.0) return dist(p, w);
    Point projection = v + t * (w - v);  // Projection falls on the segment
    return dist(p, projection);
  }

  bool line_intersection(const Point& p1, const Point& p2, const Point& p3, const Point& p4, Point* intersection) {

    double x1 = p1.x(), x2 = p2.x(), x3 = p3.x(), x4 = p4.x();
    double y1 = p1.y(), y2 = p2.y(), y3 = p3.y(), y4 = p4.y();
     
    double d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    // If d is zero, there is no intersection
    if (fzero(d)) return false;
     
    // Get the x and y
    double pre = (x1*y2 - y1*x2), post = (x3*y4 - y3*x4);
    double x = ( pre * (x3 - x4) - (x1 - x2) * post ) / d;
    double y = ( pre * (y3 - y4) - (y1 - y2) * post ) / d;
     
    // Check if the x and y coordinates are within both lines
    //if ( x < min(x1, x2) || x > max(x1, x2) ||
    //x < min(x3, x4) || x > max(x3, x4) ) return NULL;
    //if ( y < min(y1, y2) || y > max(y1, y2) ||
    //y < min(y3, y4) || y > max(y3, y4) ) return NULL;
     
    // Return the point of intersection
    if (intersection != NULL) {
      intersection->x() = x;
      intersection->y() = y;
    }
    return true;
  }

  bool segment_intersects_segment(const Point& p1, const Point& p2, const Point& p3, const Point& p4) {
    double d1 = direction(p3, p4, p1);
    double d2 = direction(p3, p4, p2);
    double d3 = direction(p1, p2, p3);
    double d4 = direction(p1, p2, p4);
  
    if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
        ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) {
        return true;
    }
      
    // Check boudary cases (end point lies on another segment)
    if (fzero(d1) && on_segment(p3, p4, p1))
        return true;
    if (fzero(d2) && on_segment(p3, p4, p2))
        return true;
    if (fzero(d3) && on_segment(p1, p2, p3))
        return true;
    if (fzero(d4) && on_segment(p1, p2, p4))
        return true;

    return false;
  }

  // return true if point a is closer to base than point b
  inline bool closer(const Point& base, const Point& a, const Point& b) {
    return (a - base).squaredNorm() < (b - base).squaredNorm();
  }

  // return true if point a is on the same side as point b with respect to the base (assume that the three points are on the same line
  inline bool same_side(const Point& base, const Point& a, const Point& b) {
    return fnneg((a.x() - base.x()) * (b.x() - base.x())) &&
           fnneg((a.y() - base.y()) * (b.y() - base.y()));
  }

  // treat the segment as an infinite line and return the updated beam
  Beam2D truncate_beam(const Beam2D& beam, const Segment& segment) {
    Beam2D out = beam;
    Point x, y;
    if (line_intersection(out.base, out.a, segment.a, segment.b, &x)) {
      if (closer(out.base, x, out.a) && same_side(out.base, x, out.a)) {
        out.a = x;
      }
    }
    if (line_intersection(out.base, out.b, segment.a, segment.b, &y)) {
      if (closer(out.base, y, out.b) && same_side(out.base, y, out.b)) {
        out.b = y;
      }
    }
    return out;
  }

  Beam2D truncate_beam(const Beam2D& beam, const vector<Segment>& segments) {
    Beam2D out = beam;
    for (int i = 0; i < segments.size(); ++i) {
      if (segment_intersects_segment(beam.base, beam.a, segments[i].a, segments[i].b) ||
          segment_intersects_segment(beam.base, beam.b, segments[i].a, segments[i].b)) {
        out = truncate_beam(out, segments[i]);
      }
    }
    return out;
  }

  void truncate_beams(const vector<Beam2D> beams, const vector<Segment>& segments, vector<Beam2D>* truncated_beams) {
    assert (truncated_beams != NULL);
    truncated_beams->clear();
    for (int i = 0; i < beams.size(); ++i) {
      truncated_beams->push_back(truncate_beam(beams[i], segments));
    }
  }

  inline double triangle_area(const Point& p0, const Point& p1, const Point& p2) {
    double area = ((p1.x() - p0.x())*(p2.y() - p0.y()) - (p2.x() - p0.x())*(p1.y()- p0.y()))/2.0;
    return fabs(area);
  }

  bool inside(const Point& p, const Beam2D& beam) {
    const Point &a = beam.a;
    const Point &b = beam.b;
    const Point &c = beam.base;

    double tot_area = triangle_area(a, b, c);
    double area1 = triangle_area(p, a, b);
    double area2 = triangle_area(p, a, c);
    double area3 = triangle_area(p, b, c);

    return feq(tot_area, area1 + area2 + area3);
  }

  bool inside(const Point& p, const vector<Beam2D>& beams) {
    assert (beams.size() > 0);
    for (int i = 0; i < beams.size(); ++i) {
      if (inside(p, beams[i])) {
        return true;
      }
    }
    return false;
  }

  // signed distance from a point to a beam
  double sgndist(const Point& p, const Beam2D& beam) {
    const Point &a = beam.a;
    const Point &b = beam.b;
    const Point &c = beam.base;
    
    double dist1 = dist_to_segment(p, a, b);
    double dist2 = dist_to_segment(p, a, c);
    double dist3 = dist_to_segment(p, b, c);
    if (inside(p, beam)) {
      return -fmin(fmin(dist1, dist2), dist3);
    } else {
      return fmin(fmin(dist1, dist2), dist3);
    }
  }

  // Note: this methods assumes that the beams are sorted counterclockwise by their angles with respect to the base
  // and that all beams share the same base point
  double sgndist(const Point& p, const vector<Beam2D>& beams) {
    assert (beams.size() > 0);
    if (!inside(p, beams)) {
      double ret = INFINITY;
      for (int i = 0; i < beams.size(); ++i) {
        ret = fmin(ret, sgndist(p, beams[i]));
      }
      return ret;
    } else {
      // only consider the distance from the point p to the enclosing segments   
      double ret = INFINITY;
      for (int i = 0; i < beams.size(); ++i) {
        ret = fmin(ret, dist_to_segment(p, beams[i].a, beams[i].b));
      }
      ret = fmin(ret, dist_to_segment(p, beams[0].base, beams[0].a));
      ret = fmin(ret, dist_to_segment(p, beams[beams.size() - 1].base, beams[beams.size() - 1].b));
      for (int i = 0; i < beams.size() - 1; ++i) {
        ret = fmin(ret, dist_to_segment(p, beams[i].b, beams[i+1].a));
      }
      return -ret;
    }
  }

  ostream& operator<< (ostream& out, const Beam2D& beam) {
    out << "Beam {base: (" << beam.base.transpose() << "), a: (" << beam.a.transpose() << "), b: (" << beam.b.transpose() << ")}";
    return out;
  }

  void partition(const Beam2D& beam, int nparts, vector<Beam2D>* beams) {
    assert (nparts > 0);
    assert (beams != NULL);
    beams->clear();
    double x1 = beam.a.x(), y1 = beam.a.y(), x2 = beam.b.x(), y2 = beam.b.y();
    double dx = (x2 - x1) / nparts, dy = (y2 - y1) / nparts;
    for (int i = 0; i < nparts; ++i) {
      Beam2D part_beam;
      part_beam.base = beam.base;
      part_beam.a = Point(x1 + dx * i, y1 + dy * i);
      part_beam.b = Point(x1 + dx * (i + 1), y1 + dy * (i + 1));
      beams->push_back(part_beam);
    }
  }

  void truncate_univariate_gaussian(double x, double cur_mean, double cur_var, double *out_mean, double *out_var) {
    double sd = sqrt(cur_var);
    double y = (x - cur_mean) / sd;
    double z = pdf(standard_normal, y) / cdf(standard_normal, y);
    *out_mean = cur_mean - z * sd;
    *out_var = cur_var * (1.0 - y*z - z*z);
  }

  void truncate_belief(const vector<Beam2D>& beams, const Vector2d& cur_mean, const Matrix2d& cur_cov, Vector2d* out_mean, Matrix2d* out_cov) {
    assert (out_mean != NULL);
    assert (out_cov != NULL);
    Point p1 = beams[0].a;
    Point p2 = beams.back().b;
    Point max_point = p1;
    double max_dist = 0.0;
    if (feq(p1.x(), p2.x())) {
      *out_mean = cur_mean;
      *out_cov = cur_cov;
      return;
    }
    for (int i = 0; i < beams.size(); ++i) {
      const Beam2D &beam = beams[i];
      double tmpdist;
      if (!segment_intersects_segment(beam.base, beam.a, p1, p2) && (tmpdist = dist_to_line(beam.a, p1, p2)) > max_dist) {
        max_dist = tmpdist;
        max_point = beam.a;
      }
      if (!segment_intersects_segment(beam.base, beam.b, p1, p2) && (tmpdist = dist_to_line(beam.b, p1, p2)) > max_dist) {
        max_dist = tmpdist;
        max_point = beam.b;
      }
    }
    double x1 = p1.x(), y1 = p1.y(), x2 = p2.x(), y2 = p2.y(), x3 = max_point.x(), y3 = max_point.y();
    Vector2d c; c << y1 - y2, x2 - x1;
    double d = (y2 - y1) * x3 - (x2 - x1) * y3;
    truncate_gaussian(c, d, cur_mean, cur_cov, out_mean, out_cov);
  }
}
