#include "TROOT.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TError.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TProfile2D.h"
#include "TChain.h"
#include "TStyle.h"
#include "TString.h"
#include "TVector3.h"
#include "TCanvas.h"
#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

#define PI 3.14159

#include "/grid/fermiapp/products/larsoft/eigen/v3_3_3/include/eigen3/Eigen/Dense"

//#include "/usr/local/Cellar/eigen/3.3.4/include/eigen3/Eigen/Dense" //Needed on MACOS
using namespace std;

struct Point {
	float x;
	float y;
	float z;
	float q;
};
struct PCAResults {
  	TVector3 centroid;
	pair<TVector3,TVector3> endPoints;
	float length;
	TVector3 eVals;
	vector<TVector3> eVecs;
};
struct TrkPoint{
    double c;
    double x;
    double y;
    double z;
    double q;
};
struct by_y { 
    bool operator()(TrkPoint const &a, TrkPoint const &b) { 
        if(a.y == b.y) return a.x > b.x;
        else return a.y > b.y;
    }
};
struct reverse_by_y { 
    bool operator()(TrkPoint const &a, TrkPoint const &b) { 
        if(a.y == b.y) return a.x < b.x;
        else return a.y < b.y;
    }
};
typedef vector<TrkPoint> track_def;
typedef vector<Point> PointCloud;
void LoadPointCloud(PointCloud &points, const track_def &ord_trk);
PCAResults DoPCA(const PointCloud &points);
double Pythagoras(double x1,double x2,double y1,double y2,double z1,double z2);
vector<double> Unit_Vec(double x1,double y1,double z1);
vector<double> Unit_Vec_NO(double x1,double x2,double y1,double y2,double z1,double z2);
double dotProdFunc(double x1,double x2,double y1,double y2,double z1,double z2);

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////MAIN PROGRAM STARTS////////////////////////////////////////

int main(int argc, char **argv){
	TFile *f_output;
	
	int track_num = atoi(argv[2]);
	//string track_id = argv[2];
	//int int_track_id = atoi(argv[2]);
	///////////////////
	//Define Parameters

	//Ordering Algorithm Parameters
	double alpha = 5.;
	double min_cosphi = -0.9;
	double phi = acos(min_cosphi) * 180. / PI;
	unsigned ang_points = 5;
	int min_points_trk = ang_points*2;
	double min_cone_ang = 0.75; //41.4 deg

	//Cluster Algorithm Parameters
	double EVR = 0.9;
	double min_costheta = -0.85;
	double max_costheta =  0.96;
	double prev_win_size = 10;
	double post_win_size = 10;

  	///////////////////
  	//Define root output file
  	f_output = TFile::Open("TEST.root","RECREATE");
  	//f_output = TFile::Open("results.root","RECREATE");
  	//f_output = TFile::Open(Form("results_alpha_%d_ang_%d.root",(int)alpha,(int)phi),"RECREATE");

  	TNtuple *nt_trk_info = new TNtuple("nt_trk_info","nt_trk_info","run_num:ev_num:cluster_id");
  	TNtuple *nt_brk_pts = new TNtuple("nt_brk_pts","nt_brk_pts","run_num:ev_num:cluster_id:alpha:dist:angle:ford_x:ford_y:ford_z:lalpha_x:lalpha_y:lalpha_z");
  	TNtuple *nt_vectors = new TNtuple("nt_vectors","nt_vectors", "run_num:ev_num:cluster_id:angle:pca_x_s:pca_y_s:pca_z_s:pca_x_e:pca_y_e:pca_z_e:lalpha_x_e:lalpha_y_e:lalpha_z_e");
  	TNtuple *nt_pca_points = new TNtuple("nt_pca_points","nt_pca_points","run_num:ev_num:cluster_id:pca_x:pca_y:pca_z");
  	


  	TNtuple *nt_trk_pts = new TNtuple("nt_trk_pts", "nt_trk_points","run_num:ev_num:cluster_id:michel:selected:ordered:x:y:z:point_ord:vertex");

  	///////////////////
  	//OUTPUT COMMENTS
  	cout << "Using only topological cuts with this selection." << endl;

  	///////////////////
  	//READ IN MICHEL LIST
  	ifstream csv_infile("Michel_candidates_vertex_v6.csv");
  	vector<string> TrackData;
  	std::vector<std::vector<double> > Michel_candidates;
  	std::string mline;
  	while (getline(csv_infile, mline,'\n')){
   		TrackData.push_back(mline); //Get each line of the file as a string
  	}
  	int s = TrackData.size();
  	for (unsigned int i=1; i<s; ++i){
   		std::vector<double> v_michel;
    	std::size_t first_comma = TrackData[i].find(",");      // position of the end of the name of each one in the respective string
    	std::size_t second_comma = TrackData[i].find(",", first_comma + 1);
    	std::size_t third_comma = TrackData[i].find(",", second_comma + 1);
    	std::size_t fourth_comma = TrackData[i].find(",", third_comma + 1);
    	std::size_t fifth_comma = TrackData[i].find(",", fourth_comma + 1);
    	double mrun = std::stod(TrackData[i].substr(0,TrackData[i].size()));
    	double meve = std::stod(TrackData[i].substr(first_comma+1,TrackData[i].size()));
    	double mtrk = std::stod(TrackData[i].substr(second_comma+1,TrackData[i].size()));
    	double mvrX = std::stod(TrackData[i].substr(third_comma+1,TrackData[i].size()));
	    double mvrY = std::stod(TrackData[i].substr(fourth_comma+1,TrackData[i].size()));
	    double mvrZ = std::stod(TrackData[i].substr(fifth_comma+1,TrackData[i].size()));
    	//if(mvrX == -1.0 && mvrY == -1.0 && mvrZ == -1.0){
    	//	continue;
	    //}
    	v_michel.push_back(mrun);
    	v_michel.push_back(meve);
    	v_michel.push_back(mtrk);
    	v_michel.push_back(mvrX);
    	v_michel.push_back(mvrY);
    	v_michel.push_back(mvrZ);
    	Michel_candidates.push_back(v_michel);
   	}
   	int mcand_size = Michel_candidates.size();
  	///////////////////////////////////////////////////////////////////////////////////////////////////
  	std::vector< std::vector<int> > kept_tracks;
  	int counter = 0;
  	int counter_x = 0;
  	int small_tracks = 0;
  	//int count = ang_count*evr_count;
	double purity, efficiency, ver_rms = 0., ver_mean = 0.;
	int total_num_tracks = 0;
	int tracks_survived_ord_alg = 0;
	int michels_survived_ord_alg = 0;
	int michel_count = 0;
	int track_selected_as_michel = 0;
	int eigenval_cut = 0;
	int ord_alg_cutouts = 0;
	int ord_alg_cutouts_michels = 0;
	std::string line;
	std::ifstream ifs(argv[1]);
		//  cout << angles[ang]  << ", " << EVRs[evr] << endl;
	while(std::getline(ifs, line)){
	  	//cout << line << endl;
	    gROOT->Reset();
	    gErrorIgnoreLevel = kError;
	    TString filename;
	    filename.Form("%s",line.c_str());    
	    TFile *infile = new TFile(filename);
	    //Extract Event Metadata
	    TTree *Trun = (TTree*)infile->Get("Trun");
	    Int_t run_num;
	    Int_t ev_num;
	    Trun->SetBranchAddress("runNo",&run_num);
	    Trun->SetBranchAddress("eventNo",&ev_num);
	    Trun->GetEntry(0);
	    //cout << "Looking at run " << run_num << " from event " << ev_num << endl;
	    //Extract Coordinate information
	    TTree *T_charge_cluster = (TTree*)infile->Get("T_charge_cluster_nfc"); 
	    Double_t cluster_id;
	    Double_t qx;
	    Double_t qy;
	    Double_t qz;
	    Double_t qc;
	    T_charge_cluster->SetBranchAddress("qx",&qx);
	    T_charge_cluster->SetBranchStatus("qx", kTRUE);
	    T_charge_cluster->SetBranchAddress("qy",&qy);
	    T_charge_cluster->SetBranchStatus("qy", kTRUE);
	    T_charge_cluster->SetBranchAddress("qz",&qz);
	    T_charge_cluster->SetBranchStatus("qz", kTRUE);
	    T_charge_cluster->SetBranchAddress("qc",&qc);
	    T_charge_cluster->SetBranchStatus("qc", kTRUE);
	    T_charge_cluster->SetBranchAddress("cluster_id", &cluster_id);
	    T_charge_cluster->SetBranchStatus("cluster_id", kTRUE);
	    int size = T_charge_cluster->GetEntries();
	    /////////////////////////////////////////////////////////////
	    //Extract Clusters///////////////////////////////////////////
	    std::vector<Int_t> clusters;
	    Int_t prev_cval;
    	for (int i = 0; i < size; ++i){
      		T_charge_cluster -> GetEntry(i);
    		if (i == 0){
        		clusters.push_back(cluster_id);
        		prev_cval = cluster_id;
      		}else if(prev_cval != cluster_id){
        		prev_cval = cluster_id;
        		clusters.push_back(cluster_id);
      		}
      	}
    	//Looking at tracks individually
    	int num_clusters = clusters.size();
   		//cout << "HERE " << num_clusters << endl;
    	//Loop through individual clusters
    	std::vector<int> event_kept_trks;
    	for (int c = 0; c < num_clusters; ++c){
    		int cluster = clusters[c];
    		cout << total_num_tracks << "\t clusters processed" << endl;
    		cout << "Looking at Run: " << run_num << ", Event: " << ev_num << ", Cluster: " << cluster << endl;
    		total_num_tracks += 1;
      		
      		string clus_id = to_string(ev_num) + to_string(cluster);
      		int int_clus_id  = stoi(clus_id);
      		//cout << clus_id << endl;

      		//PARAMETERS
      		bool is_michel = false, is_selected = false, is_ordered = false;

      		//Determine if cluster is a Michel from hand scan
      		for (int cand = 0; cand < mcand_size; ++cand){
	    			if(Michel_candidates[cand][0] == run_num && Michel_candidates[cand][1] == ev_num && Michel_candidates[cand][2] == cluster){
	    				is_michel = true;
	    				break;
	    			}
	    	}
      		nt_trk_info->Fill(run_num,ev_num,cluster);

      		track_def trk;
      		//Load every point (cluster id, x, y, z, charge) of a track into the trk object.
      		for (int i = 0; i < size; ++i){
	        	T_charge_cluster -> GetEntry(i);
	        	//Will only store information for current cluster
	        	if(cluster_id != cluster) continue;

	        	if(track_num != -1){
	        		if(int_clus_id != track_num) continue;
	        	}
		        TrkPoint tempPoint;
		        tempPoint.c = cluster_id;
		        tempPoint.x = qx;
		        tempPoint.y = qy;
		        tempPoint.z = qz;
		        tempPoint.q = qc;
		        trk.push_back(tempPoint);
	      	}

	      	//track size has to be larger than the moving window size
	      	if(trk.size() < min_points_trk + 1){
	      		small_tracks += 1;
	      		continue; // #CUT
	      	}
	      	//////////////////////////////////////////////////
	      	////////ORDERING ALGORITHM BEGINS/////////////////
	      	//Sort track in descending y value
	      	std::sort(trk.begin(), trk.end(), by_y());
	      	int trk_size = trk.size();
	      	track_def points_left;
	      	track_def points_gd;
	      	//Store points being tested by ordering algorithm
	      	for (int i = 1; i < trk_size; ++i){
	        	TrkPoint tempPoint;
	        	tempPoint.c = trk[i].c;
	        	tempPoint.x = trk[i].x;
	        	tempPoint.y = trk[i].y;
	        	tempPoint.z = trk[i].z;
	        	tempPoint.q = trk[i].q;
	        	points_left.push_back(tempPoint);
	      	}
	      	int pl_size = points_left.size();

	      	track_def ord_trk;
	      	//Highest y value point is the first point in the oredered track
	      	ord_trk.push_back(trk[0]);
		    double old_dist = 10000000.;
		    int low_dist_at = -1;
		    double dist;

		    int m = 0;
		    int i = 0;

		    double low_ord_y = 2000.;
		    double ang_eigV_x=0., ang_eigV_y=0., ang_eigV_z=0.;
		    double ob_temp_x, ob_temp_y, ob_temp_z;
		    double ftest_point_x, ftest_point_y, ftest_point_z; 
		    double pca_eVec_x, pca_eVec_y, pca_eVec_z; 
          	double ford_point_x, ford_point_y, ford_point_z, last_dist;
          	double low_old_dist = 1000.;
		    double ang_dotP, fangle;

		    std::vector<double> vec_pca_ang;
		    std::vector<double> vec_pca;
		    std::vector<double> vec_ob_temp;
		    std::vector<double> vec_temp;
		    std::vector<double> vec_last_ord;
		    double pca_end_x, pca_end_y, pca_end_z;
		    double pca_start_x, pca_start_y, pca_start_z;
		    double pca_dotP;
		    std::vector<double> pca_vec;
		    std::vector<double> cone_vec;
		    double cone_dotP;
		    double flip;


		    track_def pca_points;
		    double ang_decide;
		    TrkPoint VertexPoint;

	    	while(pl_size != 0){

	    		track_def ang_chunk;
	    		vec_pca_ang.clear();
	    		if(ord_trk.size() > ang_points){
	    			for (unsigned p = ord_trk.size() - ang_points; p <= ord_trk.size(); ++p){
		     			TrkPoint ang_tempPoint;
				        ang_tempPoint.c = ord_trk[p].c;
				        ang_tempPoint.x = ord_trk[p].x;
				        ang_tempPoint.y = ord_trk[p].y;
				        ang_tempPoint.z = ord_trk[p].z;
				        ang_tempPoint.q = ord_trk[p].q;
				        ang_chunk.push_back(ang_tempPoint);	
	    			}
		    		PointCloud ang_pointcloud;
		    		PCAResults ang_results;
		    		LoadPointCloud(ang_pointcloud, ang_chunk);
				    ang_results = DoPCA(ang_pointcloud);
					
				    vec_pca_ang.push_back(ang_results.eVecs[0](0));
					vec_pca_ang.push_back(ang_results.eVecs[0](1));
					vec_pca_ang.push_back(ang_results.eVecs[0](2));

					float xVal = ang_results.endPoints.first(0);
					float yVal = ang_results.endPoints.first(1);
  					float zVal = ang_results.endPoints.first(2);

  					int pca_counter = 0;
	  				while ((((ang_results.endPoints.first(0) < ang_results.endPoints.second(0)) && (xVal < ang_results.endPoints.second(0))) || ((ang_results.endPoints.first(0) >= ang_results.endPoints.second(0)) && (xVal > ang_results.endPoints.second(0)))) && (((ang_results.endPoints.first(1) < ang_results.endPoints.second(1)) && (yVal < ang_results.endPoints.second(1))) || ((ang_results.endPoints.first(1) >= ang_results.endPoints.second(1)) && (yVal > ang_results.endPoints.second(1)))) && (((ang_results.endPoints.first(2) < ang_results.endPoints.second(2)) && (zVal < ang_results.endPoints.second(2))) || ((ang_results.endPoints.first(2) >= ang_results.endPoints.second(2)) && (zVal > ang_results.endPoints.second(2))))) {
    					xVal += 0.5*ang_results.eVecs[0](0);
    					yVal += 0.5*ang_results.eVecs[0](1);
    					zVal += 0.5*ang_results.eVecs[0](2);
    					
    					if(pca_counter == 0){
    						pca_start_x = xVal;
    						pca_start_y = xVal;
    						pca_start_z = xVal;
    					}
    					pca_counter += 1;
  					}
  					cout << pca_counter << endl;
  					pca_end_x = xVal;
  					pca_end_y = yVal;
  					pca_end_z = zVal;
  					pca_vec = Unit_Vec_NO(pca_end_x, pca_start_x, pca_end_y, pca_start_y, pca_end_z, pca_start_z);
  					pca_dotP = dotProdFunc(pca_vec[0],vec_pca_ang[0],pca_vec[1],vec_pca_ang[1],pca_vec[2],vec_pca_ang[2]);
  					if (pca_dotP > -pca_dotP){
  						flip = 1.;
  					}else{
  						flip = -1.;
  					}

  					vec_pca_ang[0] = flip * vec_pca_ang[0];
  					vec_pca_ang[1] = flip * vec_pca_ang[1];
  					vec_pca_ang[2] = flip * vec_pca_ang[2];
			    }else{
			    	vec_pca_ang.push_back(0);
			    	vec_pca_ang.push_back(0);
			    	vec_pca_ang.push_back(0);
			    }
	    		

			    bool cone_test_fail = true;
			    bool ongoing_cone_test = true;
				bool closest_point_found = false;

			    int j = 0;
	        	//for (int j = 0; j < pl_size; ++j){
	        	while(ongoing_cone_test){
	        		TrkPoint cone_check;
	        		cone_check.x = points_left[j].x;
	        		cone_check.x = points_left[j].y;
	        		cone_check.x = points_left[j].z;

	        		cone_vec = Unit_Vec_NO(points_left[j].x,ord_trk.back().x,points_left[j].y,ord_trk.back().y,points_left[j].z,ord_trk.back().z);
	        		cone_dotP = dotProdFunc(cone_vec[0],vec_pca_ang[0],cone_vec[1],vec_pca_ang[1],cone_vec[2],vec_pca_ang[2]);

	        		if(cone_dotP > min_cone_ang){
	        			dist = Pythagoras(points_left[j].x,ord_trk.back().x,points_left[j].y,ord_trk.back().y,points_left[j].z,ord_trk.back().z);
	        			if(dist < alpha){
		          			if (dist < old_dist){
		            			old_dist = dist;
		            			low_dist_at = j;
		            			cone_test_fail = false;
		            			closest_point_found = true;
		          			}
	        			}
	        		}
	        		j += 1;
	        		if(j >= pl_size) break;
	        	}

	        	//while(cone_test_fail && ){
	        	if (cone_test_fail){
	        		for (int j = 0; j < pl_size; ++j){
		        		dist = Pythagoras(points_left[j].x,ord_trk.back().x,points_left[j].y,ord_trk.back().y,points_left[j].z,ord_trk.back().z);
		        		if(dist < alpha){
			          		if (dist < old_dist){
			            		old_dist = dist;
			            		low_dist_at = j;
			            		closest_point_found = true;
			          		}
			          	}
        			}
	        	}
	        	
	        	if(closest_point_found == false){
	        		for (int j = 0; j < pl_size; ++j){
		        		dist = Pythagoras(points_left[j].x,ord_trk.back().x,points_left[j].y,ord_trk.back().y,points_left[j].z,ord_trk.back().z);
		          		if (dist < old_dist){
		            		old_dist = dist;
		            		low_dist_at = j;
		          		}
        			}

	        	}

	        	//Point with the shortest distance to the last ordered point
	        	TrkPoint tempPoint;
	        	tempPoint.c = points_left[low_dist_at].c;
	        	tempPoint.x = points_left[low_dist_at].x;
	        	tempPoint.y = points_left[low_dist_at].y;
	        	tempPoint.z = points_left[low_dist_at].z;
	        	tempPoint.q = points_left[low_dist_at].q;

	        	ob_temp_x = tempPoint.x - ord_trk.back().x;
	        	ob_temp_y = tempPoint.y - ord_trk.back().y;
	        	ob_temp_z = tempPoint.z - ord_trk.back().z;

	        	vec_ob_temp = Unit_Vec(ob_temp_x,ob_temp_y,ob_temp_z);

	        	if(closest_point_found){
	        		//Keep the next point
	        		ord_trk.push_back(tempPoint);
	          		if (tempPoint.y < low_ord_y){
	          			low_ord_y = tempPoint.y;
	          		}
	          		old_dist = 10000000;
	          		points_left.erase(points_left.begin() + low_dist_at);
	          		pl_size = points_left.size();
	          		i = 0;
	        	}else{
	        		if(old_dist < low_old_dist){
			    		vec_last_ord.clear();
			    		vec_last_ord.push_back(ord_trk.back().x);
			    		vec_last_ord.push_back(ord_trk.back().y);
			    		vec_last_ord.push_back(ord_trk.back().z);

			    		//cout << old_dist << endl;
		                last_dist = old_dist;
		                ftest_point_x = tempPoint.x;
		                ftest_point_y = tempPoint.y;
		                ftest_point_z = tempPoint.z;

		                vec_temp = vec_ob_temp;

		                //cout << vec_temp[0] << endl;

		                pca_points = ang_chunk;
		                vec_pca.clear();
		                vec_pca.push_back(vec_pca_ang[0]);
		                vec_pca.push_back(vec_pca_ang[1]);
		                vec_pca.push_back(vec_pca_ang[2]);

		                //pca_eVec_x = ang_eigV_x;
		                //pca_eVec_y = ang_eigV_y;
		                //pca_eVec_z = ang_eigV_z;

		                //cout << acos(vec_temp[0]*ang_decide*vec_pca_ang[0] + vec_temp[1]*ang_decide*vec_pca_ang[1]) *180./PI << endl;

		                fangle = cone_dotP;
		                low_old_dist = old_dist;
              		}
			    	points_gd.push_back(tempPoint);
			        old_dist = 10000000;
			        points_left.erase (points_left.begin() + low_dist_at);
			        pl_size = points_left.size();
			        i++;
			        break;
	        	}
	        	//cout << pl_size << endl;
	        	if (pl_size == 0) break;
	    	}
			double bottom_dist;
			bottom_dist = abs(trk.back().y - low_ord_y);
			// If distance between lowest y value of unordered track and 
			// lowest y value of ordered track is greater than 10 cm.
			fangle = acos(fangle) * 180. / PI;
			//cout << vec_last_ord.size() << endl;
			if(vec_last_ord.size() != 0){
				nt_brk_pts->Fill(run_num,ev_num,cluster, alpha, last_dist, fangle, vec_last_ord[0], vec_last_ord[1], vec_last_ord[2], ftest_point_x, ftest_point_y, ftest_point_z);
				nt_vectors->Fill(run_num, ev_num, cluster, fangle, vec_last_ord[0], vec_last_ord[1], vec_last_ord[2], vec_pca[0], vec_pca[1], vec_pca[2], vec_temp[0], vec_temp[1], vec_temp[2]);
			}
			//nt_study->Fill(run_num,ev_num,cluster, alpha, last_dist, fangle, -1, -1, -1, ftest_point_x, ftest_point_y, ftest_point_z);
			for (int i = 0; i < pca_points.size(); ++i){
	        	nt_pca_points->Fill(run_num,ev_num,cluster,pca_points[i].x,pca_points[i].y,pca_points[i].z);
	        }
			if(bottom_dist > 10.){
					for (int cand = 0; cand < mcand_size; ++cand){
	    				if(Michel_candidates[cand][0] == run_num && Michel_candidates[cand][1] == ev_num && Michel_candidates[cand][2] == cluster){
	    					//cout << line << endl;
	    					ord_alg_cutouts_michels += 1;
	    					//
	    					
	    					//double norm_ord_back = Pythagoras(ord_trk.back().x,0.,ord_trk.back().y,0., ord_trk.back().z,0.);
	    					//double norm_eVec = Pythagoras(ord_trk.back().x + pca_eVec_x,0., ord_trk.back().y + pca_eVec_y,0., ord_trk.back().z + pca_eVec_x,0.);
	    					//double norm_test_point = Pythagoras(ftest_point_x,0., ftest_point_y,0., ftest_point_z,0.);
	    					//nt_study->Fill(run_num,ev_num,cluster, alpha, last_dist, fangle, ord_trk.back().x, ord_trk.back().y, ord_trk.back().z, ftest_point_x, ftest_point_y, ftest_point_z);
	    					//nt_study->Fill(run_num,ev_num,cluster, alpha, last_dist, fangle, vec_last_ord[0], vec_last_ord[1], vec_last_ord[2], ftest_point_x, ftest_point_y, ftest_point_z);
	    					//nt_vectors->Fill(run_num, ev_num, cluster, fangle, ord_trk.back().x, ord_trk.back().y, ord_trk.back().z, vec_pca[0], vec_pca[1], vec_pca[2], vec_temp[0], vec_temp[1], vec_temp[2]);
	    					//nt_vectors->Fill(run_num, ev_num, cluster, fangle, vec_last_ord[0], vec_last_ord[1], vec_last_ord[2], vec_pca[0], vec_pca[1], vec_pca[2], vec_temp[0], vec_temp[1], vec_temp[2]);
	    					//cout << run_num << ", " << ev_num << ", " << cluster << "\t Cut by Ordering Algorithm" << endl;
	    					//cout << "/Users/ivan/Work/Cosmics_Michel_Studies/Track_Plotter/Michel_Plots/Run_" << run_num << "_Event_" << ev_num << "_Cluster_" << cluster << ".pdf" << endl;
	        				//for (int i = 0; i < pca_points.size(); ++i){
	        				//	nt_pca_points->Fill(run_num,ev_num,cluster,pca_points[i].x,pca_points[i].y,pca_points[i].z);
	        				//}
	    					break;
	    				}

	    			}
				//continue; //#CUT
			}else{
				is_ordered = true;
				for (int cand = 0; cand < mcand_size; ++cand){
	    			if(Michel_candidates[cand][0] == run_num && Michel_candidates[cand][1] == ev_num && Michel_candidates[cand][2] == cluster){
	    				//cout << run_num << ", " << ev_num << ", " << cluster << "\t Cut by Ordering Algorithm" << endl;
	    				michels_survived_ord_alg += 1;
	    				break;
	    			}
	    		}
			}
			tracks_survived_ord_alg += 1;
			
			//cout << "HERE" << endl;
			//Finished Ordering Points
		    //////////////////////////
			if(is_ordered){
				//////////////////////////
				//Starting Moving Window
			 	track_def prev_chunk;
			 	track_def post_chunk;
			 	PointCloud prev_points;
			 	PointCloud post_points;
			 	//TrkPoint VertexPoint;
			 	//cout << "#################################################################" << endl;
				//cout << "Looking at track " << cluster << " from run = " << run_num << ", event = " << ev_num << endl;
				//cout << "Before window size: " << prev_win_size << "; After window size: " << post_win_size << endl;
			 	double dotProd, min_ang = -360.;
			 	double ev_lowest = 100000;
				double min_prod = 1000000000;
			 	double rr = 0.0;
				double prev_eigenratio, post_eigenratio, cross_ratios;
				double vertex_prev_ratio, vertex_post_ratio, vertex_cratio;
			 	track_def kept_prev_chunk;
			 	track_def kept_post_chunk;
			 	double prev_eVecs[3];
			 	double post_eVecs[3];
			 	int vertex;
			 	int iternum = 2400;
			 	double prev_angle;
			 	double vertex_res;
			 	unsigned trk_points = 100;
			    if(ord_trk.size() < min_points_trk*2) continue;
			    for (int i = prev_win_size + 1; i < ord_trk.size() - post_win_size - 1; ++i){
					track_def prev_chunk;
				 	track_def post_chunk;
				 	PointCloud prev_points;
				 	PointCloud post_points;    
				 	PCAResults prev_results;
				 	PCAResults post_results;
				 	TrkPoint prev_first_point;
				 	TrkPoint post_last_point;
				 	prev_first_point = ord_trk[i - prev_win_size];
				 	post_last_point = ord_trk[i + post_win_size];
					if(i == iternum) cout << "Before kink " << endl;
		     		for (int j = i - prev_win_size; j < i; ++j){
		     			TrkPoint prev_tempPoint;
				        prev_tempPoint.c = ord_trk[j].c;
				        prev_tempPoint.x = ord_trk[j].x;
				        prev_tempPoint.y = ord_trk[j].y;
				        prev_tempPoint.z = ord_trk[j].z;
				        prev_tempPoint.q = ord_trk[j].q;
				        prev_chunk.push_back(prev_tempPoint);	
						//if(i == iternum) cout << ord_trk[j].x << ", " << ord_trk[j].y << ", " << ord_trk[j].z << endl;
		     		}
					if(i == iternum) cout << "After kink" << endl;
			     	for (int j = i + 1; j < i + post_win_size + 1; ++j){
		     			TrkPoint post_tempPoint;
				        post_tempPoint.c = ord_trk[j].c;
				        post_tempPoint.x = ord_trk[j].x;
				        post_tempPoint.y = ord_trk[j].y;
				        post_tempPoint.z = ord_trk[j].z;
				        post_tempPoint.q = ord_trk[j].q;
				        post_chunk.push_back(post_tempPoint);
						//if(i == iternum) cout<< ord_trk[j].x << ", " << ord_trk[j].y << ", " << ord_trk[j].z << endl;	
			     	}
			     	LoadPointCloud(prev_points, prev_chunk);
			     	LoadPointCloud(post_points, post_chunk);
		     		prev_results = DoPCA(prev_points);
		     		post_results = DoPCA(post_points);
		     		TrkPoint prev_new_eVec_plus;
					prev_new_eVec_plus.x = ord_trk[i].x + prev_results.eVecs[0](0);
					prev_new_eVec_plus.y = ord_trk[i].y + prev_results.eVecs[0](1);
					prev_new_eVec_plus.z = ord_trk[i].z + prev_results.eVecs[0](2);
		     		TrkPoint prev_new_eVec_minus;
					prev_new_eVec_minus.x = ord_trk[i].x - prev_results.eVecs[0](0);
					prev_new_eVec_minus.y = ord_trk[i].y - prev_results.eVecs[0](1);
					prev_new_eVec_minus.z = ord_trk[i].z - prev_results.eVecs[0](2);
		     		TrkPoint post_new_eVec_plus;
					post_new_eVec_plus.x = ord_trk[i].x + post_results.eVecs[0](0);
					post_new_eVec_plus.y = ord_trk[i].y + post_results.eVecs[0](1);
					post_new_eVec_plus.z = ord_trk[i].z + post_results.eVecs[0](2);
		     		TrkPoint post_new_eVec_minus;
					post_new_eVec_minus.x = ord_trk[i].x - post_results.eVecs[0](0);
					post_new_eVec_minus.y = ord_trk[i].y - post_results.eVecs[0](1);
					post_new_eVec_minus.z = ord_trk[i].z - post_results.eVecs[0](2);
					double prev_dist_plus, prev_dist_minus, post_dist_plus, post_dist_minus;
					double prev_decide, post_decide;
					prev_dist_plus 	= Pythagoras(prev_first_point.x,prev_new_eVec_plus.x,prev_first_point.y,prev_new_eVec_plus.y,prev_first_point.z,prev_new_eVec_plus.z);
					prev_dist_minus = Pythagoras(prev_first_point.x,prev_new_eVec_minus.x, prev_first_point.y,prev_new_eVec_minus.y,prev_first_point.z,prev_new_eVec_minus.z);
					post_dist_plus 	= Pythagoras(post_last_point.x,post_new_eVec_plus.x,post_last_point.y,post_new_eVec_plus.y,post_last_point.z,post_new_eVec_plus.z);
					post_dist_minus = Pythagoras(post_last_point.x,post_new_eVec_minus.x,post_last_point.y,post_new_eVec_minus.y,post_last_point.z,post_new_eVec_minus.z);
					if(prev_dist_minus > prev_dist_plus){ 
						prev_decide = 1.;
					}else{
						prev_decide = -1.;
					}
					if(post_dist_minus > post_dist_plus){
						post_decide = 1.;
					}else{
						post_decide = -1.;
					}
					prev_results.eVecs[0] = prev_decide*prev_results.eVecs[0];
					post_results.eVecs[0] = post_decide*post_results.eVecs[0];
					prev_eigenratio = prev_results.eVals(0)/(prev_results.eVals(0) + prev_results.eVals(1) + prev_results.eVals(2));
					post_eigenratio = post_results.eVals(0)/(post_results.eVals(0) + post_results.eVals(1) + post_results.eVals(2));
		     		dotProd = prev_results.eVecs[0](0)*post_results.eVecs[0](0) + prev_results.eVecs[0](1)*post_results.eVecs[0](1) + prev_results.eVecs[0](2)*post_results.eVecs[0](2);
		     		if(i == prev_win_size){ 
		     			prev_angle = dotProd;
		     		}else{
		     			if(abs(prev_angle - dotProd) > 1.9){
		     				prev_angle = dotProd;
		     				continue;
		     			}
		     		}
					if(prev_eigenratio < EVR) continue;
					if(post_eigenratio < EVR) continue;
					if(dotProd > min_ang){
						//cout <<  acos(dotProd) * 180./PI << endl;
						min_ang = dotProd;	
						vertex_prev_ratio = prev_eigenratio;
						vertex_post_ratio = post_eigenratio;
						vertex_cratio = cross_ratios;
		     			vertex = i;
						VertexPoint.c = ord_trk[vertex].c;
		     			VertexPoint.x = ord_trk[vertex].x;
		     			VertexPoint.y = ord_trk[vertex].y;
		     			VertexPoint.z = ord_trk[vertex].z;
		     			VertexPoint.q = ord_trk[vertex].c;
		     			kept_prev_chunk = prev_chunk;
		     			kept_post_chunk = post_chunk;
		     			prev_eVecs[0] = prev_results.eVecs[0](0);
						prev_eVecs[1] = prev_results.eVecs[0](1);
						prev_eVecs[2] = prev_results.eVecs[0](2);
						post_eVecs[0] = post_results.eVecs[0](0);
						post_eVecs[1] = post_results.eVecs[0](1);
						post_eVecs[2] = post_results.eVecs[0](2);
	     			}
		    	}
		   		//////////////////////////////////////////////////////////////////////
			    double dist_low_vert_y;
			    dist_low_vert_y = abs(VertexPoint.y - low_ord_y);
			    //cout << min_ang << ", " << acos(min_ang) * 180./PI << endl;
			    //Michel electron cutoff distance based on energy spectrum
			    if(is_ordered != true) continue;
			    if(dist_low_vert_y > 15) continue;
			    if(min_ang < min_costheta) continue;
			    if(min_ang > max_costheta) continue;
			    track_selected_as_michel += 1;
			    is_selected = true;
			
		    	//cout << vertex_prev_ratio << ", " << vertex_post_ratio << ", " << min_ang << ", "<< acos(min_ang) * 180./PI << ", " << run_num << ", " << ev_num << ", " << cluster << endl;
		    	//cout << min_ang << ", " << acos(min_ang) * 180./PI << endl;
		    	for (int cand = 0; cand < mcand_size; ++cand){
		    		if(Michel_candidates[cand][0] == run_num && Michel_candidates[cand][1] == ev_num && Michel_candidates[cand][2] == cluster){
		    			//vertex_res = sqrt(pow(VertexPoint.x - Michel_candidates[cand][3], 2.) + pow(VertexPoint.y - Michel_candidates[cand][4], 2.) + pow(VertexPoint.z - Michel_candidates[cand][5], 2.));
		    			vertex_res = Pythagoras(Michel_candidates[cand][3],VertexPoint.x,Michel_candidates[cand][4],VertexPoint.y,Michel_candidates[cand][5],VertexPoint.z);
		    			ver_rms += pow(vertex_res, 2.);
		    			ver_mean += vertex_res;
		    			michel_count += 1;
		    			break;
		    		}
		    	}
		    }
		   	
	    	for(int i = 0; i < ord_trk.size(); ++i){
	    		if(ord_trk.at(i).x == VertexPoint.x && ord_trk.at(i).y == VertexPoint.y && ord_trk.at(i).z == VertexPoint.z){
	    			nt_trk_pts->Fill(run_num,ev_num,cluster,is_michel,is_selected,is_ordered,ord_trk.at(i).x,ord_trk.at(i).y,ord_trk.at(i).z,true,true);
	    		}else{
	    			nt_trk_pts->Fill(run_num,ev_num,cluster,is_michel,is_selected,is_ordered,ord_trk.at(i).x,ord_trk.at(i).y,ord_trk.at(i).z,true,false);
	    		}
	    	}
	    	
	    	for(int i = 0; i < points_gd.size(); ++i){
	    		nt_trk_pts->Fill(run_num,ev_num,cluster,is_michel,is_selected,is_ordered,points_gd.at(i).x,points_gd.at(i).y,points_gd.at(i).z,false,false);
	    	}
   		}
   		infile->Close();
	}
    purity = ((float)michel_count)/((float)track_selected_as_michel);
    efficiency = ((float)michel_count)/((float)(s-1));
    ver_rms = sqrt((1./((float)michel_count)) * ver_rms);
    ver_mean = (1./((float)michel_count)) * ver_mean;
    cout << "###################################################################" << endl;
    cout << "Total number of tracks = " << total_num_tracks << "; Number of Michels in sample = " << s - 1 << endl;
    cout << "Tracks smaller than the pca window = " << small_tracks << endl;
    cout << "###################################################################" << endl;
    cout << "################### Ordering of Tracks ############################" << endl;
    cout << "###################################################################" << endl;
    cout << "Tracks after ordering algorithm = " << tracks_survived_ord_alg << " with alpha = " << alpha << endl;
    cout << "Michel clusters that survived the ordering algorithm = " << michels_survived_ord_alg << endl;
    cout << "Michel clusters that were cut out by ordering algorithm = " << ord_alg_cutouts_michels << endl;
    cout << "Angle Criterion: Cos(/phi) < " << min_cosphi << " -> /phi > " << phi << endl;
    cout << "###################################################################" << endl;
    cout << "###################### Michel ID part #############################" << endl;
    cout << "###################################################################" << endl;
    cout << "Tracks that were selected as Michels  = " << track_selected_as_michel << endl;
    cout << "Tracks correctly selected as selected as Michels = " << michel_count << endl;
    cout << "Cuts applied on bk, ak: " << EVR << endl;
    cout << "Cuts applied on Cos(/theta): Min = " << min_costheta << ", Max = " << max_costheta << endl;
    cout << "Cuts applied on /theta: Max = " << acos(min_costheta) * 180./ PI << ", Min = " << acos(max_costheta) * 180./ PI << endl;
    cout << "###################################################################" << endl;
    cout << "########################## Reults #################################" << endl;
    cout << "###################################################################" << endl;
    cout << "Purity = " << purity << ", Efficiency = " << efficiency << ", Vertex RMS = " << ver_rms << ", Vertex Mean = " << ver_mean << endl;
    //cout << "Cuts applied on /theta: Max = " << acos(min_costheta) * 180./PI << ", Max = " << acos(max_costheta) * 180/PI << endl;
    cout << "####################################################################" << endl;


    //nt_study -> Fill(alpha,angles[ang],EVRs[evr],purity,efficiency,ver_rms,ver_mean,tracks_survived_ord_alg,track_selected_as_michel,michel_count);  
  	
  	f_output->Write();
  	f_output->Close();
  	cout << "DONE" << endl;
  	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////END OF MAIN PROGRAM/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////START OF FUNCTIONS//////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

double Pythagoras(double x1,double x2,double y1,double y2,double z1,double z2){
	double dist;
	dist = sqrt(pow(x2-x1,2.) + pow(y2-y1,2.) + pow(z2-z1,2.));
	return dist;
}

vector<double> Unit_Vec(double x1,double y1,double z1){
	std::vector<double> v;
	double norm;
	norm = Pythagoras(x1,0.0,y1,0.0,z1,0.0);
	v.push_back(x1/norm);
	v.push_back(y1/norm);
	v.push_back(z1/norm);
	return v;
}

vector<double> Unit_Vec_NO(double x1,double x2,double y1,double y2,double z1,double z2){
	std::vector<double> v;
	double norm;
	norm = Pythagoras(x1,x2,y1,y2,z1,z2);
	v.push_back((x2-x1)/norm);
	v.push_back((y2-y1)/norm);
	v.push_back((z2-z1)/norm);
	return v;
}

double dotProdFunc(double x1,double x2,double y1,double y2,double z1,double z2){
	double dotP;
	dotP = x1*x2 + y1*y2 + z1*z2;
	return dotP;
}


void LoadPointCloud(PointCloud &points, const track_def &ord_trk) {
  for (int i = 0; i < ord_trk.size(); ++i){
    Point tempPoint;
    tempPoint.x = ord_trk.at(i).x;
    tempPoint.y = ord_trk.at(i).y;
    tempPoint.z = ord_trk.at(i).z;
    tempPoint.q = ord_trk.at(i).q;
    points.push_back(tempPoint);

  }
  return;
}

PCAResults DoPCA(const PointCloud &points) {
  TVector3 outputCentroid;
  pair<TVector3,TVector3> outputEndPoints;
  float outputLength;
  TVector3 outputEigenValues;
  vector<TVector3> outputEigenVecs;
  float meanPosition[3] = {0., 0., 0.};
  unsigned int nThreeDHits = 0;
  for (unsigned int i = 0; i < points.size(); i++) {
    meanPosition[0] += points[i].x;
    meanPosition[1] += points[i].y;
    meanPosition[2] += points[i].z;
    ++nThreeDHits;
  }
  if (nThreeDHits == 0) {
    PCAResults results;
    return results; 
  }
  const float nThreeDHitsAsFloat(static_cast<float>(nThreeDHits));
  meanPosition[0] /= nThreeDHitsAsFloat;
  meanPosition[1] /= nThreeDHitsAsFloat;
  meanPosition[2] /= nThreeDHitsAsFloat;
  outputCentroid = TVector3(meanPosition[0], meanPosition[1], meanPosition[2]);
  float xi2 = 0.0;
  float xiyi = 0.0;
  float xizi = 0.0;
  float yi2 = 0.0;
  float yizi = 0.0;
  float zi2 = 0.0;
  float weightSum = 0.0;
  for (unsigned int i = 0; i < points.size(); i++) {
      const float weight(1.);
      const float x((points[i].x - meanPosition[0]) * weight);
      const float y((points[i].y - meanPosition[1]) * weight);
      const float z((points[i].z - meanPosition[2]) * weight);
      xi2  += x * x;
      xiyi += x * y;
      xizi += x * z;
      yi2  += y * y;
      yizi += y * z;
      zi2  += z * z;
      weightSum += weight * weight;
  }

  Eigen::Matrix3f sig;

  sig << xi2, xiyi, xizi,
         xiyi, yi2, yizi,
         xizi, yizi, zi2;

  sig *= 1.0 / weightSum;

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenMat(sig);

  typedef std::pair<float,size_t> EigenValColPair;
  typedef std::vector<EigenValColPair> EigenValColVector;

  EigenValColVector eigenValColVector;
  const auto &resultEigenMat(eigenMat.eigenvalues());
  eigenValColVector.emplace_back(resultEigenMat(0), 0);
  eigenValColVector.emplace_back(resultEigenMat(1), 1);
  eigenValColVector.emplace_back(resultEigenMat(2), 2);

  std::sort(eigenValColVector.begin(), eigenValColVector.end(), [](const EigenValColPair &left, const EigenValColPair &right){return left.first > right.first;} );

  outputEigenValues = TVector3(eigenValColVector.at(0).first, eigenValColVector.at(1).first, eigenValColVector.at(2).first);

  const Eigen::Matrix3f &eigenVecs(eigenMat.eigenvectors());

  for (const EigenValColPair &pair : eigenValColVector) {
     outputEigenVecs.emplace_back(eigenVecs(0, pair.second), eigenVecs(1, pair.second), eigenVecs(2, pair.second));
  }

  PCAResults results;

  Eigen::ParametrizedLine<float,3> priAxis(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)),Eigen::Vector3f(outputEigenVecs[0](0),outputEigenVecs[0](1),outputEigenVecs[0](2)));

  Eigen::Vector3f endPoint1(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)));
  Eigen::Vector3f endPoint2(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)));

  Eigen::Vector3f testPoint;
  Eigen::Vector3f projTestPoint;
  float maxDist1 = -1.0;
  float maxDist2 = -1.0;
  float dist;
  float dotP;
  for (unsigned int i = 0; i < points.size(); i++) {
    testPoint = Eigen::Vector3f(points[i].x,points[i].y,points[i].z);
    projTestPoint = priAxis.projection(testPoint);
    dist = sqrt(pow(projTestPoint(0)-outputCentroid(0),2.0)+pow(projTestPoint(1)-outputCentroid(1),2.0)+pow(projTestPoint(2)-outputCentroid(2),2.0));
    dotP = (projTestPoint(0)-outputCentroid(0))*outputEigenVecs[0](0) + (projTestPoint(1)-outputCentroid(1))*outputEigenVecs[0](1) + (projTestPoint(2)-outputCentroid(2))*outputEigenVecs[0](2);


    if ((dotP < 0.0) && (dist > maxDist1)) {
      endPoint1 = projTestPoint;
      maxDist1 = dist;
    }
    else if ((dotP > 0.0) && (dist > maxDist2)) {
      endPoint2 = projTestPoint;
      maxDist2 = dist;
    }
  }
  outputEndPoints.first = TVector3(endPoint1(0),endPoint1(1),endPoint1(2));
  outputEndPoints.second = TVector3(endPoint2(0),endPoint2(1),endPoint2(2));
  outputLength = sqrt(pow(endPoint2(0)-endPoint1(0),2.0)+pow(endPoint2(1)-endPoint1(1),2.0)+pow(endPoint2(2)-endPoint1(2),2.0));
  results.centroid = outputCentroid;
  results.endPoints = outputEndPoints;
  results.length = outputLength;
  results.eVals = outputEigenValues;
  results.eVecs = outputEigenVecs;
  return results;
}

