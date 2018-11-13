
// CollisionFree_pathplanningDlg.h: 標頭檔
//

#pragma once
#include "afxwin.h"
#include "cv.h"
#include "highgui.h"
#include "CvvImage.h"
#include "opencv2/imgproc/imgproc.hpp"
#include <vector>
using namespace cv;
using namespace std;

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
#define  total_number 5
	afx_msg void OnBnClickedButtonStart();
	afx_msg void OnBnClickedOk();
	void find_path(int x);
	void dijkstra(int source, int node_num);
	void binarization(IplImage *i_show_data, vector<vector<bool>>  &o_sca_image2);
	void find_coner(vector<vector<bool>> i_sca_image, vector <Point> &o_save_coner, int i_Interpolation);

	void trans2Voronoi(vector<vector<bool>> i_sca_image, vector <Point> i_save_coner, double(&o_Data)[8000], int i_Interpolation2);

	void Voronoi_calculate(double i_Data[8000], int x_boundary, int y_boundary, CvPoint2D64f(&o_savepoint1)[3000], CvPoint2D64f(&o_savepoint2)[3000], int &o_line_count);

	void Generalized_Voronoi(vector<vector<bool>> i_sca_image, CvPoint2D64f i_savepoint1[3000], CvPoint2D64f i_savepoint2[3000], int i_line_count, int &o_new_input_index, CvPoint2D64f(&o_new_savepoint1)[3000], CvPoint2D64f(&o_new_savepoint2)[3000]);

	void Match_point(int i_line_count, int i_new_input_index, CvPoint2D64f(&io_new_savepoint1)[3000], CvPoint2D64f(&io_new_savepoint2)[3000], float near_dis);

	void Dijkstra_path_planning(int i_highest_robot_priority, CvPoint2D64f i_robot_start[total_number], CvPoint2D64f  i_robot_end[total_number], CvPoint2D64f i_new_savepoint1[3000], CvPoint2D64f i_new_savepoint2[3000], int i_new_input_index, vector <CPoint> &o_all_point_map, vector <CvPoint2D64f> &o_all_point_map_original);

	void Path_Optimization(vector<vector<bool>> i_sca_image, vector <CvPoint2D64f> i_all_point_map_original, int i_robot_num, vector <int> &o_path_optimization);

	void MultiRobot_Path_simulation(CDC* i_pDC, IplImage * i_draw_data, vector <CPoint> i_host_path, vector<vector<bool>>  i_sca_image, vector <Point> i_save_coner, CvPoint2D64f i_robot_start_point[total_number], CvPoint2D64f i_robot_end_point[total_number], int i_car_density, vector <CPoint> &o_sim_path, vector <draw_car> &o_sim_car, IplImage *&offline_show);
	//輸入路徑、車子間距，輸出模擬路徑、車子

	void ServantRobot_Path_simulation(vector <CPoint> i_Servant_path, int i_robot_num, CvPoint2D64f &o_ServantRobot_pos, CvPoint2D64f i_robot_start_point, double &io_zdir);

	void servant_path(int i_robot_num, IplImage * i_pGrayImg, vector<vector<bool>>  i_sca_image, vector <Point> i_save_coner, CvPoint2D64f i_robot_start_point[total_number], CvPoint2D64f i_robot_end_point[total_number], vector <CPoint> i_sim_path, int i_carsize, CvVideoWriter *i_slave_recoder, int i_arrive_goal_times[total_number]);

	void Control_Methods(int control_type, double i_rho, double i_alpha, double i_beta, double i_phi, double &o_vr, double &o_vl, int &o_state);

	void simulation_car(IplImage * &live_show, CvPoint2D64f i_robot_start_point[total_number], double i_robot_zdir[total_number], int i_car_num, int i_carsize, int i_arrive_goal_times[total_number]);

	void car_crash_detector(CvPoint2D64f i_robot_start_point[total_number], int i_carsize, int &io_crash_times, char io_crash_num[total_number][total_number]);

	void car_front_detector(CvPoint2D64f i_robot_start_point[total_number], double i_robot_zdir[total_number], double i_detect_angle, int i_detect_radius, CvPoint(&o_detect_area)[total_number][2]);

	void car_path_efficiency(CvPoint2D64f i_robot_start_point, CvPoint2D64f i_robot_end_point, int reentry_tag, int old_reentry_tag, int i_car_num);

	CStatic m_show;
	int m_coner_count;
	int m_total_time;
	CStatic m_show2;
};



class float_point
{
	float x;
	float y;
};

class CElement : public CObject
{
	DECLARE_SERIAL(CElement)

protected:
	COLORREF m_Color;                       // Color of an element
	CRect m_EnclosingRect;                  // Rectangle enclosing an element
	int m_Pen;                              // Pen width

public:
	virtual ~CElement() {}                   // Virtual destructor

											 // Virtual draw operation
	virtual void Draw(CDC* pDC, const CElement* pElement = 0) const {}
	virtual void Move(const CSize& Size) {} // Move an element
	CRect GetBoundRect() const;             // Get the bounding rectangle for an element

	virtual void Serialize(CArchive& ar);      // Serialize function for CElement

protected:
	CElement() {}                            // Default constructor
};


// Class defining a vor_point object
class CVor_Point : public CElement
{
	DECLARE_SERIAL(CVor_Point)

public:
	int GetPositionY();
	int GetPositionX();
	// Function to display a line
	virtual void Draw(CDC* pDC, const CElement* pElement = 0) const;
	virtual void Move(const CSize& aSize);       // Function to move an element

												 // Constructor for a line object
	CVor_Point(const CPoint& Start, const COLORREF& Color, const int& PenWidth);

	virtual void Serialize(CArchive& ar);      // Serialize function for CVor_Point

protected:
	CPoint m_StartPoint;          // Start point of line

	CVor_Point() {}             // Default constructor - should not be used
};

class SK_point
{
public:
	int x;
	int y;
	float distant;
};

