
// CollisionFree_pathplanningDlg.h: 標頭檔
//

#pragma once
#include "afxwin.h"
#include "cv.h"
#include "highgui.h"
#include "CvvImage.h"
#include "opencv2/imgproc/imgproc.hpp"
#include <cstdio>
#include <vector>

#include "GraphElements.h"
#include "Graph.h"
#include "YenTopKShortestPathsAlg.h"
#include "DijkstraShortestPathAlg.h"

#include <boost/polygon/voronoi.hpp>
#include <boost/polygon/polygon.hpp>
// using boost::polygon::voronoi_builder;
// using boost::polygon::voronoi_diagram;
// using boost::polygon::x;
// using boost::polygon::y;
// using boost::polygon::low;
// using boost::polygon::high;

using namespace boost::polygon;
using namespace cv;
using namespace std;

typedef double coordinate_type;
typedef point_data<coordinate_type> point_type;
typedef segment_data<coordinate_type> segment_type;
typedef rectangle_data<coordinate_type> rect_type;
typedef voronoi_diagram<coordinate_type> VD;
typedef VD::cell_type cell_type;
typedef VD::cell_type::source_index_type source_index_type;
typedef VD::cell_type::source_category_type source_category_type;
typedef VD::edge_type edge_type;


class CCollisionFreepathplanningDlgAutoProxy;


struct draw_car
{
	CvPoint car[6];
	vector <CvPoint>sim_path;
};

// CCollisionFreepathplanningDlg 對話方塊
class CCollisionFreepathplanningDlg : public CDialogEx
{
	DECLARE_DYNAMIC(CCollisionFreepathplanningDlg);
	friend class CCollisionFreepathplanningDlgAutoProxy;

// 建構
public:
	CCollisionFreepathplanningDlg(CWnd* pParent = nullptr);	// 標準建構函式
	virtual ~CCollisionFreepathplanningDlg();

// 對話方塊資料
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_COLLISIONFREE_PATHPLANNING_DIALOG };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV 支援


// 程式碼實作
protected:
	CCollisionFreepathplanningDlgAutoProxy* m_pAutoProxy;
	HICON m_hIcon;

	BOOL CanExit();

	// 產生的訊息對應函式
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	afx_msg void OnClose();
	virtual void OnOK();
	virtual void OnCancel();
	DECLARE_MESSAGE_MAP()
public:
#define  total_number 5  //設定有幾台機器人+1
#define  total_target 9  //設定有多少個目標點+1
	afx_msg void OnBnClickedButtonStart();
	afx_msg void OnBnClickedOk();
	void binarization(IplImage *i_show_data, vector<vector<bool>>  &o_sca_image2);
	void find_coner(vector<vector<bool>> i_sca_image, vector <Point> &o_save_coner, int i_Interpolation);

	void trans2Voronoi(vector<vector<bool>> i_sca_image, vector <Point> i_save_coner, double(&o_Data)[3000], int i_Interpolation2);

	void boost_Voronoi_calculate(int x_boundary, int y_boundary, CvPoint2D64f(&o_savepoint1)[1500], CvPoint2D64f(&o_savepoint2)[1500], int &o_line_count);

	void clip_infinite_edge(const edge_type& edge, std::vector<point_type>* clipped_edge);

	void iterate_primary_edges1(const voronoi_diagram<double>& vd, CvPoint2D64f(&o_savepoint1)[1500], CvPoint2D64f(&o_savepoint2)[1500], int & o_line_count);

	void Generalized_Voronoi(vector<vector<bool>> i_sca_image, CvPoint2D64f i_savepoint1[1500], CvPoint2D64f i_savepoint2[1500], int i_line_count, int &o_new_input_index, CvPoint2D64f(&o_new_savepoint1)[1500], CvPoint2D64f(&o_new_savepoint2)[1500]);

	void Match_point(int i_line_count, int i_new_input_index, CvPoint2D64f(&io_new_savepoint1)[1500], CvPoint2D64f(&io_new_savepoint2)[1500], float near_dis);

	void Dijkstra_path_planning(int i_highest_robot_priority, CvPoint2D64f i_robot_start[total_number], CvPoint2D64f  i_robot_end[total_number], CvPoint2D64f i_new_savepoint1[1500], CvPoint2D64f i_new_savepoint2[1500], int i_new_input_index, vector <CPoint> &o_all_point_map, vector <CvPoint2D64f> &o_all_point_map_original);

	void Path_Optimization_20181127(vector<vector<bool>> i_sca_image, vector <CvPoint2D64f> i_all_point_map_original, int i_robot_num, vector <int> &o_path_optimization);

	void Check_Thepath_Collision(vector<vector<bool>> i_sca_image, vector <CvPoint2D64f> i_all_point_map_original, int i_robot_num, vector <CPoint> i_all_point_map, vector <int> i_path_optimization[total_number]);

	void MultiRobot_Path_simulation(IplImage * i_draw_data, vector <CPoint> i_host_path, vector<vector<bool>>  i_sca_image, vector <Point> i_save_coner, CvPoint2D64f i_robot_start_point[total_number], CvPoint2D64f i_robot_end_point[total_number], int i_car_density, vector <CPoint> &o_sim_path, vector <draw_car> &o_sim_car, IplImage *&offline_show);
	//輸入路徑、車子間距，輸出模擬路徑、車子

	void ServantRobot_Path_simulation(vector <CPoint> i_Servant_path, int i_robot_num, CvPoint2D64f &o_ServantRobot_pos, CvPoint2D64f i_robot_start_point, double &io_zdir);

	void Servant_Path(int i_robot_num, IplImage * i_pGrayImg, vector<vector<bool>>  i_sca_image, vector <Point> i_save_coner, CvPoint2D64f i_robot_start_point[total_number], CvPoint2D64f i_robot_end_point[total_number], vector <CPoint> i_sim_path, int i_carsize, CvVideoWriter *i_slave_recoder, int i_arrive_goal_times[total_number]);

	void Control_Methods(int control_type, double i_rho, double i_alpha, double i_beta, double i_phi, double &o_vr, double &o_vl, int &o_state);

	void simulation_car(IplImage * &live_show, CvPoint2D64f i_robot_start_point[total_number], double i_robot_zdir[total_number], int i_car_num, int i_carsize, int i_arrive_goal_times[total_number]);

	void car_crash_detector(CvPoint2D64f i_robot_start_point[total_number], int i_carsize, int &io_crash_times, char io_crash_num[total_number][total_number]);

	void car_front_detector(CvPoint2D64f i_robot_start_point[total_number], double i_robot_zdir[total_number], double i_detect_angle, int i_detect_radius, int i_carsize, CvPoint(&o_detect_area)[total_number][2]);

	void car_path_efficiency(CvPoint2D64f i_robot_start_point, CvPoint2D64f i_robot_end_point, int reentry_tag, int old_reentry_tag, int i_car_num);

	float my_cross(CPoint o, CPoint a, CPoint b)
	{
		return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
	}

	float my_cross2(CPoint a, CPoint b)
	{
		return (a.x) * (b.y) - (a.y) * (b.x);
	}

	bool my_intersect1D(float a1, float a2, float b1, float b2)
	{
		if (a1 > a2)
			swap(a1, a2);
		if (b1 > b2)
			swap(b1, b2);
		return max(a1, b1) <= min(a2, b2);
	}

	bool is_intersect(CPoint a1, CPoint a2, CPoint b1, CPoint b2)
	{
		return my_intersect1D(a1.x, a2.x, b1.x, b2.x)
			&& my_intersect1D(a1.y, a2.y, b1.y, b2.y)
			&& my_cross(a1, a2, b1) * my_cross(a1, a2, b2) <= 0
			&& my_cross(b1, b2, a1) * my_cross(b1, b2, a2) <= 0;
	}

	CPoint intersection_point(CPoint a1, CPoint a2, CPoint b1, CPoint b2)
	{
		CPoint a = a2 - a1, b = b2 - b1, s = b1 - a1;
		a1.x = a1.x + a.x * my_cross2(s, b) / my_cross2(a, b);
		a1.y = a1.y + a.y * my_cross2(s, b) / my_cross2(a, b);
		// 計算交點
		return a1;
	}

	CStatic m_show;
	int m_coner_count;
	int m_total_time;
	CStatic m_show2;
	BOOL check_path_change;
	CButton check_path_change2;
};


