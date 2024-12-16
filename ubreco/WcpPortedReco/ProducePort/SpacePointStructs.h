struct TrecSpacePoint {
  double x;
  double y;
  double z;
  double q;
  double nq;
  int cluster_id;
  int real_cluster_id;
  int sub_cluster_id;
};

struct TrecchargeSpacePoint {
  double x;
  double y;
  double z;
  double q;
  double nq;
  int cluster_id;
  int real_cluster_id;
  int sub_cluster_id;
  double chi2;
  double ndf;
  double pu;
  double pv;
  double pw;
  double pt;
  double reduced_chi2;
  int flag_vertex;
  int flag_shower;
  double rr;
};

struct TrecchargeblobSpacePoint {
  double x;
  double y;
  double z;
  double q;
  double nq;
  int cluster_id;
  int real_cluster_id;
  int sub_cluster_id;
  double chi2;
  double ndf;
};

struct TclusterSpacePoint {
  double x;
  double y;
  double z;
  double q;
  double nq;
  int cluster_id;
  int real_cluster_id;
  int sub_cluster_id;
  int time_slice;
  int ch_u;
  int ch_v;
  int ch_w;
};