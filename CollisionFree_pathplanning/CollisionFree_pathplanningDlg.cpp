
// CollisionFree_pathplanningDlg.cpp: 實作檔案
//

#include "stdafx.h"
#include "CollisionFree_pathplanning.h"
#include "CollisionFree_pathplanningDlg.h"
#include "DlgProxy.h"
#include "afxdialogex.h"
#include "fstream"
#include <windows.h>
#include "Voronoi.cpp"
#include "math.h"
#include "queue"
#include "CvvImage.cpp"
#include "omp.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// 對 App About 使用 CAboutDlg 對話方塊

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

	// 對話方塊資料
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_ABOUTBOX };
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支援

// 程式碼實作
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(IDD_ABOUTBOX)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// CCollisionFreepathplanningDlg 對話方塊


IMPLEMENT_DYNAMIC(CCollisionFreepathplanningDlg, CDialogEx);

CCollisionFreepathplanningDlg::CCollisionFreepathplanningDlg(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_COLLISIONFREE_PATHPLANNING_DIALOG, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
	m_pAutoProxy = nullptr;
}

CCollisionFreepathplanningDlg::~CCollisionFreepathplanningDlg()
{
	// 如果有此對話方塊的 Automation Proxy，
	//  請將其指向此對話方塊的返回指標設為 null，讓其知道
	// 所以會知道是否已經刪除對話方塊。
	if (m_pAutoProxy != nullptr)
		m_pAutoProxy->m_pDialog = nullptr;
}

void CCollisionFreepathplanningDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_STATIC_show, m_show);
	DDX_Control(pDX, IDC_STATIC_show2, m_show2);
}

BEGIN_MESSAGE_MAP(CCollisionFreepathplanningDlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_CLOSE()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_BUTTON_Start, &CCollisionFreepathplanningDlg::OnBnClickedButtonStart)
	ON_BN_CLICKED(IDOK, &CCollisionFreepathplanningDlg::OnBnClickedOk)
END_MESSAGE_MAP()


// CCollisionFreepathplanningDlg 訊息處理常式

BOOL CCollisionFreepathplanningDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// 將 [關於...] 功能表加入系統功能表。

	// IDM_ABOUTBOX 必須在系統命令範圍之中。
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != nullptr)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// 設定此對話方塊的圖示。當應用程式的主視窗不是對話方塊時，
	// 框架會自動從事此作業
	SetIcon(m_hIcon, TRUE);			// 設定大圖示
	SetIcon(m_hIcon, FALSE);		// 設定小圖示

	// TODO: 在此加入額外的初始設定

	return TRUE;  // 傳回 TRUE，除非您對控制項設定焦點
}

void CCollisionFreepathplanningDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

// 如果將最小化按鈕加入您的對話方塊，您需要下列的程式碼，
// 以便繪製圖示。對於使用文件/檢視模式的 MFC 應用程式，
// 框架會自動完成此作業。

void CCollisionFreepathplanningDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // 繪製的裝置內容

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// 將圖示置中於用戶端矩形
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// 描繪圖示
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

// 當使用者拖曳最小化視窗時，
// 系統呼叫這個功能取得游標顯示。
HCURSOR CCollisionFreepathplanningDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

// 如果控制器仍保持其物件之一，而使用者
// 關閉 UI 時，Automation 伺服器不應該結束。
// 這些訊息處理常式會確認是否仍在使用 Proxy，
// 如果已結束使用，則會隱藏 UI，但是對話方塊
// 仍保持在附近。

int car_crash_times = 0;
int d[2000];       // 紀錄起點到各個點的最短路徑長度
int parent[2000];  // 紀錄各個點在最短路徑樹上的父親是誰
int w[2000][2000];    // 一張有權重的圖
bool visit[2000];  // 紀錄各個點是不是已在最短路徑樹之中
double robot_zdir[total_number] = { 0 };
double show_path_efficiency[total_number] = { 0 };
vector <int> show_path;
vector <CvPoint> all_robot_path[total_number];
vector <CvPoint2D64f> all_robot_path_double[total_number];
vector <double> all_robot_dir[total_number];
CvPoint2D64f fix_robot_start_point[total_number] = { 0 }, fix_robot_end_point[total_number] = { 0 };
CvPoint2D64f delay_own_obstacle[total_number][2] = { 0 };
bool stop_program = false;
bool car_front_num[total_number][total_number] = { false };
char crash_num[total_number][total_number] = { 0 };
voronoi_diagram<double> vd;
rect_type brect_;
vector<point_type> point_data_;
vector<segment_type> segment_data_;

struct voro_Point {
	int a;
	int b;
	voro_Point(int x, int y) : a(x), b(y) {}
};

struct Segment {
	voro_Point p0;
	voro_Point p1;
	Segment(int x1, int y1, int x2, int y2) : p0(x1, y1), p1(x2, y2) {}
};

namespace boost {
	namespace polygon {

		template <>
		struct geometry_concept<voro_Point> {
			typedef point_concept type;
		};

		template <>
		struct point_traits<voro_Point>
		{
			typedef int coordinate_type;

			static inline coordinate_type get(const voro_Point& point, orientation_2d orient)
			{
				return (orient == HORIZONTAL) ? point.a : point.b;
			}
		};

		template <>
		struct geometry_concept<Segment> {
			typedef segment_concept type;
		};

		template <>
		struct segment_traits<Segment> {
			typedef int coordinate_type;
			typedef voro_Point point_type;

			static inline point_type get(const Segment& segment, direction_1d dir) {
				return dir.to_int() ? segment.p1 : segment.p0;
			}
		};
	}  // polygon
}  // boost

segment_type retrieve_segment(const cell_type& cell)
{
	source_index_type index = cell.source_index() - point_data_.size();
	return segment_data_[index];
}

point_type retrieve_point(const cell_type& cell)
{
	source_index_type index = cell.source_index();
	source_category_type category = cell.source_category();
	if (category == SOURCE_CATEGORY_SINGLE_POINT)
	{
		return point_data_[index];
	}
	index -= point_data_.size();
	if (category == SOURCE_CATEGORY_SEGMENT_START_POINT)
	{
		return low(segment_data_[index]);
	}
	else
	{
		return high(segment_data_[index]);
	}
}

void CCollisionFreepathplanningDlg::OnBnClickedButtonStart()
{

	CWnd* CW_vo = (CWnd *)GetDlgItem(IDC_STATIC_show2);
	CDC* pDC2 = CW_vo->GetWindowDC();

	remove("路徑輸出.txt");
	remove("子機器人輸出.txt");
	remove("master_recoder.avi");
	remove("slave_recoder_1.avi");
	remove("slave_recoder_2.avi");
	remove("slave_recoder_3.avi");
	remove("slave_recoder_4.avi");


	IplConvKernel *pKernel_small = NULL;
	IplConvKernel *pKernel_small2 = NULL;
	pKernel_small = cvCreateStructuringElementEx(11, 11, 6, 6, CV_SHAPE_RECT, NULL);
	vector<vector<bool>>  sca_image;  //縮圖後二值化的結果
	vector <int> path_optimization;
	vector <Point> save_coner;
	vector <CPoint> all_point_map;
	vector <CPoint>  jump_path_optimization, jump_path_optimization_simulation;
	vector <CvPoint2D64f> all_point_map_original;
	vector <draw_car> car_simulation;
	int corner_count = 0;
	int line_count = 0;  //有多少VD線段
	int new_input_index = 0;
	int path_optimization_size_change = 0;
	bool image_change = 0;
	double Data[8000];
	//	template <double T, size_t N>
	CElement* pElement = 0;
	CvPoint2D32f *point2 = 0;
	CvPoint2D64f savepoint1[3000] = { 0.0 }, savepoint2[3000] = { 0.0 };
	CvPoint2D64f new_savepoint1[3000] = { 0.0 }, new_savepoint2[3000] = { 0.0 };
	LARGE_INTEGER tStart, tEnd, ts;

	IplImage * pGrayImg = NULL;
	IplImage * draw_data = NULL;
	IplImage * check_change = NULL;
	IplImage * show_data = NULL; //縮小十倍的矩陣
	IplImage * read_data_old = NULL;
	IplImage * read_data = NULL;

	read_data = cvLoadImage("格點化地圖.bmp", 0);
	read_data_old = cvCreateImage(cvGetSize(read_data), read_data->depth, 1);
	pGrayImg = cvCreateImage(cvGetSize(read_data), read_data->depth, 1);
	draw_data = cvCreateImage(cvGetSize(read_data), read_data->depth, 3);
	check_change = cvCreateImage(cvGetSize(read_data), read_data->depth, 1);
	show_data = cvCreateImage(cvSize(88, 88), IPL_DEPTH_8U, 1);
	m_show2.SetWindowPos(&wndTop, 10, 10, draw_data->width, draw_data->height, SWP_SHOWWINDOW);
	char path0[100];
	int photo_conunt = 0;

	CvPoint2D64f robot_start_point[total_number] = { 0 }, robot_end_point[total_number] = { 0 };  //宣告多機器人之起點與終點

// 	robot_start_point[0].x = 180;  //路徑起始與終點，請參照圖片給定
// 	robot_start_point[0].y = 820;
// 	robot_end_point[0].x = 640;
// 	robot_end_point[0].y = 39;

// 	robot_start_point[1].x = 150;  //路徑起始與終點，請參照圖片給定
// 	robot_start_point[1].y = 150;
// 	robot_end_point[1].x = 700;
// 	robot_end_point[1].y = 700;

	robot_start_point[1].x = 50;  //路徑起始與終點，請參照圖片給定
	robot_start_point[1].y = 350;
	robot_end_point[1].x = 800;
	robot_end_point[1].y = 800;
	if (total_number > 2)
	{
		// 		robot_start_point[2].x = 700;  //路徑起始與終點，請參照圖片給定
		// 		robot_start_point[2].y = 700;
		// 		robot_end_point[2].x = 150;
		// 		robot_end_point[2].y = 150;

		robot_start_point[2].x = 800;  //路徑起始與終點，請參照圖片給定
		robot_start_point[2].y = 530;
		robot_end_point[2].x = 50;
		robot_end_point[2].y = 135;
	}

	if (total_number > 3)
	{
		// 		robot_start_point[3].x = 150;  //路徑起始與終點，請參照圖片給定
		// 		robot_start_point[3].y = 700;
		// 		robot_end_point[3].x = 700;
		// 		robot_end_point[3].y = 150;

		robot_start_point[3].x = 250;  //路徑起始與終點，請參照圖片給定
		robot_start_point[3].y = 800;
		robot_end_point[3].x = 550;
		robot_end_point[3].y = 50;
	}

	if (total_number > 4)
	{
		// 		robot_start_point[4].x = 700;  //路徑起始與終點，請參照圖片給定
		// 		robot_start_point[4].y = 150;
		// 		robot_end_point[4].x = 150;
		// 		robot_end_point[4].y = 700;

		robot_start_point[4].x = 800;  //路徑起始與終點，請參照圖片給定
		robot_start_point[4].y = 250;
		robot_end_point[4].x = 100;
		robot_end_point[4].y = 650;
	}
	memcpy(fix_robot_start_point, robot_start_point, sizeof(fix_robot_start_point));
	memcpy(fix_robot_end_point, robot_end_point, sizeof(fix_robot_end_point));

	photo_conunt = 1;
	while (true)
	{
		QueryPerformanceFrequency(&ts);
		QueryPerformanceCounter(&tStart);

		sprintf_s(path0, "photo\\L%d.png", photo_conunt);
		fstream in_image0(path0, ios::in);

		if (!in_image0)
			break;

		read_data = cvLoadImage(path0, 0);

		if (1)
		{
			//-------------------------清除數據------------------------------
			save_coner.clear();
			sca_image.clear();
			path_optimization.clear();
			all_point_map.clear();
			all_point_map_original.clear();
			show_path.clear();
			memset((unsigned char*)draw_data->imageData, 0, draw_data->imageSize);
			memset(d, 0, sizeof(d));
			memset(parent, 0, sizeof(parent));
			memset(visit, 0, sizeof(visit));
			memset(w, 0, sizeof(w));

			//---------------------------------------------------------------------

			cvResize(read_data, pGrayImg, CV_INTER_NN);//讀黑白影像用的
	//		cvCvtColor(read_data, pGrayImg, CV_RGB2GRAY);  //讀彩色影像用的
			cvErode(pGrayImg, pGrayImg, pKernel_small, 3);  //侵蝕的相反(因為是白底)
			cvDilate(pGrayImg, pGrayImg, pKernel_small, 2);  //膨脹的相反
			cvCvtColor(pGrayImg, draw_data, CV_GRAY2RGB);
			//			cvSaveImage("給連通物件用的.bmp", pGrayImg);

						//數值要依據縮小倍率與格點pixel數決定


			cvResize(pGrayImg, show_data, CV_INTER_NN);


			//輸入圖片，輸出二值資料
			binarization(show_data, sca_image);
			//輸入二值資料，輸出角點
			find_coner(sca_image, save_coner, 4);
			//將角點轉換為準備要丟入Voronoi運算的格式
			trans2Voronoi(sca_image, save_coner, Data, 8);
			//計算狹義Voronoi，輸入角點資料與邊界，輸出兩個矩陣 

// 			Voronoi_calculate(Data, show_data->width, show_data->height, savepoint1, savepoint2, line_count);  //20180926拿掉Master
// 			//計算廣義Voronoi，待改
// 			Generalized_Voronoi(sca_image, savepoint1, savepoint2, line_count, new_input_index, new_savepoint1, new_savepoint2);
// 			//VD點會破碎，將其重新聚合
// 			Match_point(line_count, new_input_index, new_savepoint1, new_savepoint2, 2);
// 			//Dijkstra路徑搜尋，輸入點連接資訊跟數量
// 			Dijkstra_path_planning(0, robot_start_point, robot_end_point, new_savepoint1, new_savepoint2, new_input_index, all_point_map, all_point_map_original);
// 			//路徑優化，輸入二值資訊與原本路徑
// 			Path_Optimization(sca_image, all_point_map_original, path_optimization);



			//-------------------------------------------繪圖---------------------------------------   //20180926拿掉Master
//			cvDilate(draw_data, draw_data, pKernel_small, 2);  //膨脹的相反 
// 			for (int i = 0; i < save_coner.size(); i++)  //角點圖
// 			{
// 				cvLine(draw_data, cvPoint(save_coner[i].x * 10, save_coner[i].y * 10), cvPoint(save_coner[i].x * 10, save_coner[i].y * 10), CV_RGB(0, 250, 250), 8);
// 			}
// 			for (int i = 0; i < line_count; i++)   //VD圖
// 			{
// 				cvLine(draw_data, cvPoint(savepoint1[i].x * 10, savepoint1[i].y * 10), cvPoint(savepoint2[i].x * 10, savepoint2[i].y * 10), CV_RGB(0, 0, 255), 1);
// 			}

// 			for (int i = 0; i < new_input_index; i++)  //GVD圖
// 			{
// 				cvLine(draw_data, cvPoint(new_savepoint1[i].x * 10, new_savepoint1[i].y * 10), cvPoint(new_savepoint2[i].x * 10, new_savepoint2[i].y * 10), CV_RGB(250, 100, 100), 1);
// 			}
// 
// 			for (int path_index = 0; path_index < show_path.size() - 1; path_index++) //畫出路徑圖
// 			{
// 				cvLine(draw_data, cvPoint(all_point_map[show_path[path_index]].x, all_point_map[show_path[path_index]].y), cvPoint(all_point_map[show_path[path_index + 1]].x, all_point_map[show_path[path_index + 1]].y), CV_RGB(0, 0, 255), 3);
// 			}
// 
// 			float pathpoint_dis = 200;
// 			int jump_num2 = 20000, jump_num;
// 			int path_opt = 0;
// 
// 			for (int path_opt = 0; path_opt < path_optimization.size(); path_opt++) //路徑優化搜尋
// 			{
// 				if (path_optimization_size_change != path_optimization.size())
// 				{
// 					jump_num = 3000000;
// 				}
// 
// 				if (pathpoint_dis > 16)  //16
// 				{
// 					//					first_in = sqrtf(pow((all_point_map[path_optimization[0]].x - all_point_map[path_optimization[1]].x), 2) + pow((all_point_map[path_optimization[0]].y - all_point_map[path_optimization[1]].y), 2)) / 2;
// 					jump_num2 = path_optimization[path_opt];
// 
// 					if (jump_num2 != jump_num)
// 						jump_path_optimization.push_back(all_point_map[path_optimization[path_opt]]);
// 
// 				}
// 				else
// 				{
// 					jump_num = path_optimization[1];
// 					//					jump_path_optimization_save = all_point_map[path_optimization[1]];
// 
// 				}
// 				if (path_opt == path_optimization.size() - 1)
// 					break;
// 				else
// 				{
// 					pathpoint_dis = sqrtf(pow((all_point_map[path_optimization[0]].x - all_point_map[path_optimization[path_opt + 1]].x), 2) + pow((all_point_map[path_optimization[0]].y - all_point_map[path_optimization[path_opt + 1]].y), 2));
// 					cvLine(draw_data, cvPoint(all_point_map[path_optimization[path_opt]].x, all_point_map[path_optimization[path_opt]].y), cvPoint(all_point_map[path_optimization[path_opt + 1]].x, all_point_map[path_optimization[path_opt + 1]].y), CV_RGB(0, 150, 0), 3);
// 					//					cvCircle(draw_data, cvPoint(all_point_map[path_optimization[path_opt]].x, all_point_map[path_optimization[path_opt]].y), 7, CV_RGB(0, 150, 0), 2);
// 				}
// 			}
// 
// 			path_optimization_size_change = path_optimization.size();
// 			//			jump_path_optimization.push_back(all_point_map[path_optimization[path_opt]]);
// 
// 
// 
// 
// 			cvSaveImage("大圖輸出.png", draw_data);
			//-------------------------------------------繪製小地圖---------------------------------------

// 			IplImage * itest2 = NULL;
// 			IplImage * itest = NULL;
// 			itest2 = cvCreateImage(cvSize(88, 88), IPL_DEPTH_8U, 1);
// 			itest = cvCreateImage(cvSize(88, 88), IPL_DEPTH_8U, 3);
// 			cvResize(read_data, itest2, CV_INTER_NN);
// 			cvCvtColor(itest2, itest, CV_GRAY2RGB);
// 
// 			for (int i = 0; i < save_coner.size(); i++)  //角點圖
// 			{
// 				cvLine(itest, cvPoint(save_coner[i].x , save_coner[i].y ), cvPoint(save_coner[i].x, save_coner[i].y), CV_RGB(0, 250, 250), 1);
// 			}
// 			for (int i = 0; i < line_count; i++)   //VD圖
// 			{
// 				cvLine(itest, cvPoint(savepoint1[i].x, savepoint1[i].y), cvPoint(savepoint2[i].x, savepoint2[i].y), CV_RGB(250, 200, 100), 1);
// 			}
// 			for (int i = 0; i < new_input_index; i++)  //GVD圖
// 			{
// 				cvLine(itest, cvPoint(new_savepoint1[i].x, new_savepoint1[i].y), cvPoint(new_savepoint2[i].x, new_savepoint2[i].y), CV_RGB(250, 100, 100), 1);
// 			}
// 
// 			for (int path_index = 0; path_index < show_path.size() - 1; path_index++) //畫出路徑圖
// 			{
// 				cvLine(itest, cvPoint(all_point_map[show_path[path_index]].x, all_point_map[show_path[path_index]].y), cvPoint(all_point_map[show_path[path_index + 1]].x, all_point_map[show_path[path_index + 1]].y), CV_RGB(0, 0, 255), 1);
// 			}
// 
// 			for (int path_opt = 0; path_opt < path_optimization.size() - 1; path_opt++) //畫出路徑優化圖
// 			{
// 
// 				cvLine(itest, cvPoint(all_point_map[path_optimization[path_opt]].x, all_point_map[path_optimization[path_opt]].y), cvPoint(all_point_map[path_optimization[path_opt + 1]].x, all_point_map[path_optimization[path_opt + 1]].y), CV_RGB(0, 150, 0), 1);
// 
// 				cvCircle(itest, cvPoint(all_point_map[path_optimization[path_opt]].x, all_point_map[path_optimization[path_opt]].y), 7, CV_RGB(0, 150, 0), 1);
// 
// 			}
// 
// 			cvSaveImage("輸出.png", itest);
// 			cvReleaseImage(&itest);
// 			cvReleaseImage(&itest2);
			//-------------------------------------------繪圖---------------------------------------

			MultiRobot_Path_simulation(pDC2, show_data, jump_path_optimization, sca_image, save_coner, robot_start_point, robot_end_point, 2, jump_path_optimization_simulation, car_simulation, draw_data);  //開始模擬
			jump_path_optimization.clear();

			QueryPerformanceCounter(&tEnd);

			m_total_time = 1000 / ((tEnd.QuadPart - tStart.QuadPart) * 1000 / (double)(ts.QuadPart));
			m_coner_count = Data[0] + 1; //顯示角點數量
			UpdateData(FALSE);
		}



		memset((unsigned char*)show_data->imageData, 0, show_data->imageSize);
		cvReleaseImage(&read_data);


		break;
		//		photo_conunt++;
	}

	OnBnClickedOk();
}


void CCollisionFreepathplanningDlg::binarization(IplImage * i_show_data, vector<vector<bool>> &o_sca_image2)
{
	vector<bool> sca_image1;

	for (int y = 0; y < i_show_data->height - 1; y++)
	{
		for (int x = 0; x < i_show_data->width - 1; x++)
		{
			int temp_of_image = cvGetReal2D(i_show_data, y, x);
			if (temp_of_image == 255)
				sca_image1.push_back(0);
			else
				sca_image1.push_back(1);
		}
		o_sca_image2.push_back(sca_image1);
		sca_image1.clear();
	}

	for (int x = 0; x < i_show_data->width - 1; x++)
	{
		sca_image1.push_back(0);
	}

	o_sca_image2.push_back(sca_image1);
	sca_image1.clear();

}

void CCollisionFreepathplanningDlg::find_coner(vector<vector<bool>> i_sca_image, vector <Point> &o_save_coner, int i_Interpolation)
{
	int result_out;
	int my_mask[3][3] = { { 0, -1, 0 },{ -1, 4, -1 },{ 0, -1, 0 } };
	int Interpolation = 0;

	for (unsigned int y = 1; y < i_sca_image.size() - 1; y++)
	{
		for (unsigned int x = 1; x < i_sca_image[0].size() - 1; x++)
		{

			result_out =
				my_mask[0][0] * i_sca_image[y - 1][x - 1] +
				my_mask[0][1] * i_sca_image[y - 1][x] +
				my_mask[0][2] * i_sca_image[y - 1][x + 1] +
				my_mask[1][0] * i_sca_image[y][x - 1] +
				my_mask[1][1] * i_sca_image[y][x] +
				my_mask[1][2] * i_sca_image[y][x + 1] +
				my_mask[2][0] * i_sca_image[y + 1][x - 1] +
				my_mask[2][1] * i_sca_image[y + 1][x] +
				my_mask[2][2] * i_sca_image[y + 1][x + 1];


			if (result_out == 1)
			{
				Interpolation++;
				if (Interpolation == i_Interpolation)
				{
					o_save_coner.push_back(cvPoint(x, y));
					Interpolation = 0;
				}
			}
			else
				Interpolation = 0;


			if (my_mask[0][1] * i_sca_image[y - 1][x] + my_mask[2][1] * i_sca_image[y + 1][x] == -2)  //濾掉直線
				continue;
			if (my_mask[1][0] * i_sca_image[y][x - 1] + my_mask[1][2] * i_sca_image[y][x + 1] == -2)  //濾掉橫線
				continue;

			if (result_out == 2 || result_out == 3)
				o_save_coner.push_back(cvPoint(x, y));
		}
	}


	for (unsigned int x = 1; x < i_sca_image[0].size() - 1; x++)
	{
		for (unsigned int y = 1; y < i_sca_image.size() - 1; y++)
		{

			result_out =
				my_mask[0][0] * i_sca_image[y - 1][x - 1] +
				my_mask[0][1] * i_sca_image[y - 1][x] +
				my_mask[0][2] * i_sca_image[y - 1][x + 1] +
				my_mask[1][0] * i_sca_image[y][x - 1] +
				my_mask[1][1] * i_sca_image[y][x] +
				my_mask[1][2] * i_sca_image[y][x + 1] +
				my_mask[2][0] * i_sca_image[y + 1][x - 1] +
				my_mask[2][1] * i_sca_image[y + 1][x] +
				my_mask[2][2] * i_sca_image[y + 1][x + 1];


			if (result_out == 1)
			{
				Interpolation++;
				if (Interpolation == i_Interpolation)
				{
					o_save_coner.push_back(cvPoint(x, y));
					Interpolation = 0;
				}
			}
			else
			{
				Interpolation = 0;
			}

		}
	}

}

void CCollisionFreepathplanningDlg::trans2Voronoi(vector<vector<bool>> i_sca_image, vector<Point> i_save_coner, double(&o_Data)[8000], int i_Interpolation2)
{
	unsigned int input_Data = 0;
	int jump_count = 0;
	jump_count = (i_sca_image.size() + 1) / i_Interpolation2;
	

	o_Data[0] = i_save_coner.size();
	for (input_Data = 1; input_Data < i_save_coner.size() + 1; input_Data++)
	{
		//		cvLine(RGB_show_data, cvPoint(save_coner[input_Data - 1].x, save_coner[input_Data - 1].y),cvPoint(save_coner[input_Data - 1].x, save_coner[input_Data - 1].y), CV_RGB(0, 255, 0), 1);
		o_Data[2 * input_Data - 1] = i_save_coner[input_Data - 1].x;
		o_Data[2 * input_Data] = i_save_coner[input_Data - 1].y;

		point_data_.push_back(voro_Point(i_save_coner[input_Data - 1].x, i_save_coner[input_Data - 1].y));
	}



	for (unsigned int i = 0; i <= i_sca_image.size(); i = i + i_Interpolation2)
	{
		for (unsigned int j = 0; j <= i_sca_image[0].size() + 1; j = j + i_Interpolation2)
		{
			if (i == 0 || i == i_sca_image.size() || j == 0 || j == i_sca_image[0].size() + 1)
			{
				point_data_.push_back(voro_Point(i, j));
				o_Data[2 * input_Data - 1] = i;
				o_Data[2 * input_Data] = j;
				input_Data++;
				o_Data[0]++;
				//				cvLine(RGB_show_data, cvPoint(i, j),cvPoint(i, j), CV_RGB(200, 200, 0), 1);
			}
		}
	}
}

void CCollisionFreepathplanningDlg::Voronoi_calculate(double i_Data[8000], int x_boundary, int y_boundary, CvPoint2D64f(&o_savepoint1)[3000], CvPoint2D64f(&o_savepoint2)[3000], int &o_line_count)
{
	//	remove("原始VD座標輸出.txt");
	//	fstream app_VD_output2("原始VD座標輸出.txt", ios::app);

	CWnd* CW_vo = (CWnd *)GetDlgItem(IDC_STATIC_show);
	CDC* pDC = CW_vo->GetWindowDC();
	CVoronoi* vor;
	vor = new CVoronoi();
	POSITION aPos;
	Site** EdgeSite;

	double *x1, *x2, *y1, *y2, let00 = 0, pxmax = x_boundary - 1, pymax = -y_boundary + 1, temp_y1, temp_x1, temp_x2, temp_y2, pymax_invers = y_boundary - 1;
	int i_line_count = 0, line_count2 = 0;
	int pos;

	// 	CPen aPen;
	// 	aPen.CreatePen(PS_SOLID, 2, RGB(255, 0, 0));
	// 	pDC->SelectObject(&aPen);
	vor->SetPoints(i_Data);
	vor->DrawEdges(pDC, 0, 1);
	CList<Edge, Edge&> *VorEdgeList = vor->GetEdges();
	CList<Edge, Edge&> *VorLineList = vor->GetLines();
	if (!VorEdgeList->IsEmpty())
	{
		aPos = VorEdgeList->GetHeadPosition();
		do {
			Edge& aEdge = VorEdgeList->GetNext(aPos);
			pos = aEdge.edgenbr;
			x1 = &aEdge.ep[0]->coord.x;
			y1 = &aEdge.ep[0]->coord.y;
			x2 = &aEdge.ep[1]->coord.x;
			y2 = &aEdge.ep[1]->coord.y;

			if (aEdge.ep[0] && aEdge.ep[1])
			{
				//				app_VD_output2 << *x1 << ", " << *y1 << " 到 " << *x2 << ", " << *y2 << ",第 " << line_count2 << endl;
				line_count2++;
			}

			if (!aEdge.ep[0] || *x1 < let00 || *y1 < let00)
			{
				if (aEdge.b != 0)
				{
					x1 = &let00;
					temp_y1 = (aEdge.c - aEdge.a*(*x1)) / aEdge.b;
					y1 = &temp_y1;
				}
				else
				{
					temp_x1 = aEdge.c / aEdge.a;
					x1 = &temp_x1;
					y1 = &let00;
				}
			}

			if (!aEdge.ep[1] || *y2 < let00 || *x2 < let00)
			{
				if (aEdge.b != 0)
				{
					x2 = &pxmax;
					temp_y2 = (aEdge.c - aEdge.a*(*x2)) / aEdge.b;
					y2 = &temp_y2;
				}
				else
				{
					temp_x2 = aEdge.c / aEdge.a;
					x2 = &temp_x2;
					y2 = &pymax;
				}
			}

			if (!*y1)
			{
				//				app_VD_output << (int)*x1 << ", " << (int)*y1 << " 到 " << (int)*x2 << ", " << (int)*y2 << "--------y1沒值---------" << endl;
				y1 = &pymax_invers;
			}

			if ((int)*x1 == (int)*x2 && (int)*y1 == (int)*y2)
				continue;

			if (*y2 < 0)
				*y2 = 0;

			// 			pDC->MoveTo((int)*x1, (int)*y1);
			// 			pDC->LineTo((int)*x2, (int)*y2);


			o_savepoint1[i_line_count].x = (*x1);
			o_savepoint1[i_line_count].y = (*y1);
			o_savepoint2[i_line_count].x = (*x2);
			o_savepoint2[i_line_count].y = (*y2);
			i_line_count++;

		} while (aPos);

	}

	o_line_count = i_line_count;

	//	app_VD_output2.close();
	ReleaseDC(pDC);
	delete vor;

}

void CCollisionFreepathplanningDlg::boost_Voronoi_calculate(int x_boundary, int y_boundary, CvPoint2D64f(&o_savepoint1)[3000], CvPoint2D64f(&o_savepoint2)[3000], int & o_line_count)
{
	o_line_count = 0;
	construct_voronoi(point_data_.begin(), point_data_.end(), &vd);
	iterate_primary_edges1(vd, o_savepoint1, o_savepoint2, o_line_count);

}


void CCollisionFreepathplanningDlg::clip_infinite_edge(const edge_type& edge, std::vector<point_type>* clipped_edge)
{
	const cell_type& cell1 = *edge.cell();
	const cell_type& cell2 = *edge.twin()->cell();
	point_type origin, direction;
	// Infinite edges could not be created by two segment sites.
	if (cell1.contains_point() && cell2.contains_point())
	{
		point_type p1 = retrieve_point(cell1);
		point_type p2 = retrieve_point(cell2);
		origin.x((p1.x() + p2.x()) * 0.5);
		origin.y((p1.y() + p2.y()) * 0.5);
		direction.x(p1.y() - p2.y());
		direction.y(p2.x() - p1.x());
	}
	else
	{
		origin = cell1.contains_segment() ?
			retrieve_point(cell2) :
			retrieve_point(cell1);
		segment_type segment = cell1.contains_segment() ? retrieve_segment(cell1) : retrieve_segment(cell2);
		coordinate_type dx = high(segment).x() - low(segment).x();
		coordinate_type dy = high(segment).y() - low(segment).y();
		if ((low(segment) == origin) ^ cell1.contains_point())
		{
			direction.x(dy);
			direction.y(-dx);
		}
		else
		{
			direction.x(-dy);
			direction.y(dx);
		}
	}
	coordinate_type side = xh(brect_) - xl(brect_);
	coordinate_type koef = side / (std::max)(fabs(direction.x()), fabs(direction.y()));
	if (edge.vertex0() == NULL)
	{
		clipped_edge->push_back(point_type(
			origin.x() - direction.x() * koef,
			origin.y() - direction.y() * koef));
	}
	else
	{
		clipped_edge->push_back(point_type(edge.vertex0()->x(), edge.vertex0()->y()));
	}
	if (edge.vertex1() == NULL)
	{
		clipped_edge->push_back(point_type(
			origin.x() + direction.x() * koef,
			origin.y() + direction.y() * koef));
	}
	else
	{
		clipped_edge->push_back(point_type(edge.vertex1()->x(), edge.vertex1()->y()));
	}
}


// Traversing Voronoi edges using edge iterator.
void CCollisionFreepathplanningDlg::iterate_primary_edges1(const voronoi_diagram<double>& vd, CvPoint2D64f(&o_savepoint1)[3000], CvPoint2D64f(&o_savepoint2)[3000], int & o_line_count)
{
	int result = 0;

	for (voronoi_diagram<double>::const_edge_iterator it = vd.edges().begin(); it != vd.edges().end(); ++it)
	{
		if (it->is_primary())
		{
			++result;

			if (result % 2 == 0)
				continue;

			std::vector<point_type> samples;
			if (!it->is_finite())
			{
				clip_infinite_edge(*it, &samples);
//				printf("infinite (%0.2lf %0.2lf) to (%0.2lf %0.2lf)\n", samples[0].x(), samples[0].y(), samples[1].x(), samples[1].y());
				o_savepoint1[o_line_count].x = samples[0].x();
				o_savepoint1[o_line_count].y = samples[0].y();
				o_savepoint2[o_line_count].x = samples[1].x();
				o_savepoint2[o_line_count].y = samples[1].y();
				o_line_count++;
			}
			else
			{
//				printf("(%0.2lf %0.2lf) to (%0.2lf %0.2lf)\n", it->vertex0()->x(), it->vertex0()->y(), it->vertex1()->x(), it->vertex1()->y());
				o_savepoint1[o_line_count].x = it->vertex0()->x();
				o_savepoint1[o_line_count].y = it->vertex0()->y();
				o_savepoint2[o_line_count].x = it->vertex1()->x();
				o_savepoint2[o_line_count].y = it->vertex1()->y();
				o_line_count++;
			}
		}

	}
}

void CCollisionFreepathplanningDlg::Generalized_Voronoi(vector<vector<bool>> i_sca_image, CvPoint2D64f i_savepoint1[3000], CvPoint2D64f i_savepoint2[3000], int i_line_count, int &o_new_input_index, CvPoint2D64f(&o_new_savepoint1)[3000], CvPoint2D64f(&o_new_savepoint2)[3000])
{
	//	remove("廣義VD座標輸出.txt");
	//	fstream app_VD_output("廣義VD座標輸出.txt", ios::app);

	int cheak_point = 0, line_distant, new_input_index = 0;
	CvPoint2D64f cutout_point[2000] = { 0 };

	for (cheak_point = 0; cheak_point < i_line_count; cheak_point++)
	{
		if (i_savepoint1[cheak_point].x == 0 && i_savepoint1[cheak_point].y == 0 && i_savepoint2[cheak_point].x == 0 && i_savepoint2[cheak_point].y == 0)
			continue;

		if (i_savepoint1[cheak_point].x >= 87)  //拯救右邊界
			i_savepoint1[cheak_point].x--;
		if (i_savepoint2[cheak_point].x >= 87)
			i_savepoint2[cheak_point].x--;
		if (i_savepoint1[cheak_point].y <= 1)
			i_savepoint1[cheak_point].y += 2;
		if (i_savepoint2[cheak_point].y <= 1)
			i_savepoint2[cheak_point].y += 2;

		line_distant = sqrt(pow(i_savepoint1[cheak_point].x - i_savepoint2[cheak_point].x, 2) + pow(i_savepoint1[cheak_point].y - i_savepoint2[cheak_point].y, 2));
		//		line_distant++;

		int cutout = 0;
		for (cutout = 0; cutout < line_distant + 2; cutout++)
		{
			//計算斜線上的點
			if (i_savepoint1[cheak_point].x >= i_savepoint2[cheak_point].x && i_savepoint1[cheak_point].y >= i_savepoint2[cheak_point].y)
			{
				cutout_point[cutout].x = i_savepoint2[cheak_point].x + (sqrt(pow((i_savepoint1[cheak_point].x - i_savepoint2[cheak_point].x), 2)) / (double)line_distant)*(double)cutout;
				cutout_point[cutout].y = i_savepoint2[cheak_point].y + (sqrt(pow((i_savepoint1[cheak_point].y - i_savepoint2[cheak_point].y), 2)) / (double)line_distant)*(double)cutout;
			}

			if (i_savepoint2[cheak_point].x >= i_savepoint1[cheak_point].x && i_savepoint1[cheak_point].y >= i_savepoint2[cheak_point].y)
			{
				cutout_point[cutout].x = i_savepoint1[cheak_point].x + (sqrt(pow((i_savepoint1[cheak_point].x - i_savepoint2[cheak_point].x), 2)) / (double)line_distant)*(double)cutout;
				cutout_point[cutout].y = i_savepoint1[cheak_point].y - (sqrt(pow((i_savepoint1[cheak_point].y - i_savepoint2[cheak_point].y), 2)) / (double)line_distant)*(double)cutout;
			}

			if (i_savepoint1[cheak_point].x >= i_savepoint2[cheak_point].x && i_savepoint2[cheak_point].y >= i_savepoint1[cheak_point].y)
			{
				cutout_point[cutout].x = i_savepoint2[cheak_point].x + (sqrt(pow((i_savepoint1[cheak_point].x - i_savepoint2[cheak_point].x), 2)) / (double)line_distant)*(double)cutout;
				cutout_point[cutout].y = i_savepoint2[cheak_point].y - (sqrt(pow((i_savepoint1[cheak_point].y - i_savepoint2[cheak_point].y), 2)) / (double)line_distant)*(double)cutout;
			}

			if (i_savepoint2[cheak_point].x >= i_savepoint1[cheak_point].x && i_savepoint2[cheak_point].y >= i_savepoint1[cheak_point].y)
			{
				cutout_point[cutout].x = i_savepoint1[cheak_point].x + (sqrt(pow((i_savepoint1[cheak_point].x - i_savepoint2[cheak_point].x), 2)) / (double)line_distant)*(double)cutout;
				cutout_point[cutout].y = i_savepoint1[cheak_point].y + (sqrt(pow((i_savepoint1[cheak_point].y - i_savepoint2[cheak_point].y), 2)) / (double)line_distant)*(double)cutout;
			}
		}

		if (cutout < 3)
		{
			// 			if ((cutout_point[0].y > i_sca_image.size() - 1) || (cutout_point[0].x > i_sca_image[0].size() - 1) || (cutout_point[0].y < 1 || (cutout_point[0].x < 1)))
			// 				continue;

// 			if (i_sca_image[round(cutout_point[0].y)][round(cutout_point[0].x)] == 0 ||/**/
// 				i_sca_image[(int)cutout_point[0].y + 1][(int)cutout_point[0].x] == 0 ||
// 				i_sca_image[(int)cutout_point[0].y - 1][(int)cutout_point[0].x] == 0 ||
// 				i_sca_image[(int)cutout_point[0].y][(int)cutout_point[0].x - 1] == 0 ||
// 				i_sca_image[(int)cutout_point[0].y][(int)cutout_point[0].x + 1] == 0/**/)
// 			{
// 				o_new_savepoint1[new_input_index] = i_savepoint1[cheak_point];
// 				o_new_savepoint2[new_input_index] = i_savepoint2[cheak_point];
// 				new_input_index++;
// 			}

		}

		for (int i = 0; i < cutout - 1; i++)
		{
			if ((cutout_point[i].y > i_sca_image.size() - 1) || (cutout_point[i].x > i_sca_image[0].size() - 1) || (cutout_point[i].y < 1 || (cutout_point[i].x < 1)))
				continue;

			if (i_sca_image[round(cutout_point[i].y)][round(cutout_point[i].x)] == 1 ||
				_isnan(cutout_point[i].y)
				/*||
				i_sca_image[round(cutout_point[i].y + 1)][round(cutout_point[i].x)] == 1 ||
				i_sca_image[round(cutout_point[i].y - 1)][round(cutout_point[i].x)] == 1||
				i_sca_image[round(cutout_point[i].y)][round(cutout_point[i].x - 1)] == 1 ||
				 i_sca_image[round(cutout_point[i].y)][round(cutout_point[i].x + 1)] == 1*/)
			{

				memset(i_savepoint1, 0, sizeof(i_savepoint1));
				memset(i_savepoint2, 0, sizeof(i_savepoint2));
				break;
			}

			if (i == cutout - 2)
			{
				o_new_savepoint1[new_input_index] = i_savepoint1[cheak_point];
				o_new_savepoint2[new_input_index] = i_savepoint2[cheak_point];

				//				app_VD_output << o_new_savepoint1[i].x * (double)10.0 << ", " << o_new_savepoint1[i].y * (double)10.0 << " 到 " << o_new_savepoint2[i].x * (double)10.0 << ", " << o_new_savepoint2[i].y * (double)10.0 << ", 第 " << new_input_index << endl;
				new_input_index++;  //廣義VD線段數量+1
				memset(i_savepoint1, 0, sizeof(i_savepoint1));
				memset(i_savepoint2, 0, sizeof(i_savepoint2));
				break;
			}
		}
	}
	o_new_input_index = new_input_index;  //輸出最後有多少廣義VD線段
}

void CCollisionFreepathplanningDlg::Match_point(int i_line_count, int i_new_input_index, CvPoint2D64f(&io_new_savepoint1)[3000], CvPoint2D64f(&io_new_savepoint2)[3000], float near_dis)
{
	int small_loop1, small_loop2, small_index = 0;
	double small_dis1, small_dis2, small_dis3;
	CvPoint2D64f temp_point;
	vector <CvPoint2D64f> center_point, center_point2;

	for (small_loop1 = 0; small_loop1 < i_line_count; small_loop1++)
	{
		for (small_loop2 = 0; small_loop2 < i_line_count; small_loop2++)
		{
			if (io_new_savepoint1[small_loop1].x == io_new_savepoint2[small_loop2].x&&
				io_new_savepoint1[small_loop1].y == io_new_savepoint2[small_loop2].y)
				continue;


			small_dis1 = sqrt(pow(io_new_savepoint1[small_loop1].x - io_new_savepoint2[small_loop2].x, 2) + pow(io_new_savepoint1[small_loop1].y - io_new_savepoint2[small_loop2].y, 2));
			if (small_dis1 < near_dis)
			{
				temp_point.x = (io_new_savepoint1[small_loop1].x + io_new_savepoint2[small_loop2].x) / 2;
				temp_point.y = (io_new_savepoint1[small_loop1].y + io_new_savepoint2[small_loop2].y) / 2;
				//io_new_savepoint2[small_loop2] = io_new_savepoint1[small_loop1];
				center_point.push_back(cvPoint2D64f(temp_point.x, temp_point.y));
			}
		}
	}


	for (small_loop1 = 0; small_loop1 < i_line_count; small_loop1++)
	{
		for (small_loop2 = 0; small_loop2 < i_line_count; small_loop2++)
		{
			small_dis1 = sqrt(pow(io_new_savepoint1[small_loop1].x - io_new_savepoint1[small_loop2].x, 2) + pow(io_new_savepoint1[small_loop1].y - io_new_savepoint1[small_loop2].y, 2));
			small_dis2 = sqrt(pow(io_new_savepoint2[small_loop1].x - io_new_savepoint2[small_loop2].x, 2) + pow(io_new_savepoint2[small_loop1].y - io_new_savepoint2[small_loop2].y, 2));
			if (small_dis1 < near_dis)
			{
				temp_point.x = (io_new_savepoint1[small_loop1].x + io_new_savepoint1[small_loop2].x) / 2;
				temp_point.y = (io_new_savepoint1[small_loop1].y + io_new_savepoint1[small_loop2].y) / 2;
			}
			if (small_dis2 < near_dis)
			{
				temp_point.x = (io_new_savepoint2[small_loop1].x + io_new_savepoint2[small_loop2].x) / 2;
				temp_point.y = (io_new_savepoint2[small_loop1].y + io_new_savepoint2[small_loop2].y) / 2;
			}
		}
	}

	center_point2.push_back(center_point[0]);

	for (unsigned int i = 1; i < center_point.size(); i++)
	{

		for (unsigned int j = 0; j < center_point2.size(); j++)
		{
			small_dis3 = sqrt(pow(center_point[i].x - center_point2[j].x, 2) + pow(center_point[i].y - center_point2[j].y, 2));

			if (i == j)
				continue;

			if (small_dis3 < near_dis)
			{
				//				center_point2.push_back(center_point[i]);
				small_index = 0;
				break;
			}
			small_index = 1;
		}

		if (small_index == 1)
		{
			small_index = 0;
			center_point2.push_back(center_point[i]);
		}

	}

	for (int i = 0; i < i_new_input_index; i++)
	{

		for (unsigned int j = 0; j < center_point2.size(); j++)
		{
			small_dis1 = sqrt(pow(io_new_savepoint1[i].x - center_point2[j].x, 2) + pow(io_new_savepoint1[i].y - center_point2[j].y, 2));
			small_dis2 = sqrt(pow(io_new_savepoint2[i].x - center_point2[j].x, 2) + pow(io_new_savepoint2[i].y - center_point2[j].y, 2));

			if (small_dis1 < near_dis)
			{
				io_new_savepoint1[i] = center_point2[j];
			}

			if (small_dis2 < near_dis)
			{
				io_new_savepoint2[i] = center_point2[j];
			}
		}
	}
}

void CCollisionFreepathplanningDlg::Dijkstra_path_planning(int i_robot_num, CvPoint2D64f i_robot_start[total_number], CvPoint2D64f  i_robot_end[total_number], CvPoint2D64f i_new_savepoint1[3000], CvPoint2D64f i_new_savepoint2[3000], int i_new_input_index, vector <CPoint> &o_all_point_map, vector <CvPoint2D64f> &o_all_point_map_original)
{
	vector <CPoint> CPoint_savepoint1;
	vector <CPoint> CPoint_savepoint2;

	for (int i = 0; i < i_new_input_index; i++)
	{
		CPoint_savepoint1.push_back(CPoint(round(i_new_savepoint1[i].x * 10), round(i_new_savepoint1[i].y * 10)));
		CPoint_savepoint2.push_back(CPoint(round(i_new_savepoint2[i].x * 10), round(i_new_savepoint2[i].y * 10)));
	}

	int loop1, loop2;
	bool onestime;
	int put_index = 0;

	//先丟入第一組，以免檢查時vector時沒東西
	o_all_point_map.push_back(CPoint_savepoint1[0]);
	o_all_point_map_original.push_back(i_new_savepoint1[0]);
	/*	app_path_num << o_all_point_map[put_index].x << ", " << o_all_point_map[put_index].y << " 第 " << put_index << endl;*/
	put_index++;

	for (int i = 0; i < i_new_input_index; i++)  //過濾重複點
	{
		onestime = 1;
		for (unsigned int j = 0; j < o_all_point_map.size(); j++)
		{
			if (o_all_point_map[j] == CPoint_savepoint1[i])
				onestime = 0;
		}
		if (onestime == 1)
		{
			o_all_point_map.push_back(CPoint_savepoint1[i]);
			o_all_point_map_original.push_back(i_new_savepoint1[i]); // 原始的也排列一次
//			app_path_num << o_all_point_map[put_index].x << ", " << o_all_point_map[put_index].y << " 第 " << put_index << endl;
			put_index++;
		}
	}

	for (int i = 0; i < i_new_input_index; i++)
	{
		onestime = 1;
		for (unsigned int j = 0; j < o_all_point_map.size(); j++)
		{
			if (o_all_point_map[j] == CPoint_savepoint2[i])
				onestime = 0;
		}
		if (onestime == 1)
		{
			o_all_point_map.push_back(CPoint_savepoint2[i]);
			o_all_point_map_original.push_back(i_new_savepoint2[i]); // 原始的也排列一次
//			app_path_num << o_all_point_map[put_index].x << ", " << o_all_point_map[put_index].y << " 第 " << put_index << endl;
			put_index++;
		}
	}

	int savetemp_index1, savetemp_index2;
	int the_point_index[2];
	int serch_point[4];
	serch_point[2] = 100000;
	serch_point[3] = 100000;

	for (int serch = 0; serch < put_index; serch++)  //此刻位置與目標點皆要搜尋最靠近的GVD點做連接
	{

		serch_point[0] = sqrt(pow((o_all_point_map[serch].x - i_robot_start[i_robot_num].x), 2) + pow((o_all_point_map[serch].y - i_robot_start[i_robot_num].y), 2));
		serch_point[1] = sqrt(pow((o_all_point_map[serch].x - i_robot_end[i_robot_num].x), 2) + pow((o_all_point_map[serch].y - i_robot_end[i_robot_num].y), 2));


		if (serch_point[2] > serch_point[0])
		{
			serch_point[2] = serch_point[0];
			the_point_index[0] = serch;
		}
		if (serch_point[3] > serch_point[1])
		{
			serch_point[3] = serch_point[1];
			the_point_index[1] = serch;
		}
	}

	o_all_point_map.push_back(CPoint(i_robot_start[i_robot_num].x, i_robot_start[i_robot_num].y));
	CPoint_savepoint1.push_back(CPoint(i_robot_start[i_robot_num].x, i_robot_start[i_robot_num].y));
	o_all_point_map_original.push_back(cvPoint2D64f(i_robot_start[i_robot_num].x / 10, i_robot_start[i_robot_num].y / 10));
	CPoint_savepoint2.push_back(CPoint(o_all_point_map[the_point_index[0]].x, o_all_point_map[the_point_index[0]].y));
	i_new_input_index++;
	put_index++;

	for (loop1 = 0; loop1 < i_new_input_index; loop1++)
	{
		for (loop2 = 0; loop2 < put_index; loop2++)
		{

			if (CPoint_savepoint1[loop1] == o_all_point_map[loop2])
				savetemp_index1 = loop2;

			if (CPoint_savepoint2[loop1] == o_all_point_map[loop2])
				savetemp_index2 = loop2;

			if (w[loop1][loop2] == 0)
				w[loop1][loop2] = 50000;

		}

		w[savetemp_index1][savetemp_index2] = sqrt(pow(CPoint_savepoint1[loop1].x - CPoint_savepoint2[loop1].x, 2) + pow(CPoint_savepoint1[loop1].y - CPoint_savepoint2[loop1].y, 2));
		w[savetemp_index2][savetemp_index1] = w[savetemp_index1][savetemp_index2];
	}


	dijkstra(put_index - 1, o_all_point_map.size());
	find_path(the_point_index[1]);

}

void CCollisionFreepathplanningDlg::Path_Optimization_20181127(vector<vector<bool>> i_sca_image, vector<CvPoint2D64f> i_all_point_map_original, int i_robot_num, vector<int>& o_path_optimization)
{

	int path_inside, break_index, line_distant, crash_index;
	unsigned int cheak_x, cheak_y, cheak_collision_index = 0, final_noncollision_index = 0;
	CvPoint2D64f cutout_point[500];
	CvPoint2D64f temp;
	vector <bool> allpath_crash_index;

	temp.x = fix_robot_end_point[i_robot_num].x / 10;
	temp.y = fix_robot_end_point[i_robot_num].y / 10;

	o_path_optimization.push_back(show_path[0]);
	show_path.push_back(i_all_point_map_original.size());
	i_all_point_map_original.push_back(temp);

	for (unsigned int i = 0; i < show_path.size() - 1;)
	{
		crash_index = 0;
		for (unsigned int j = i + 1; j < show_path.size(); j++)
		{

			line_distant = sqrt(pow(i_all_point_map_original[show_path[i]].x - i_all_point_map_original[show_path[j]].x, 2) + pow(i_all_point_map_original[show_path[i]].y - i_all_point_map_original[show_path[j]].y, 2));

			int cutout = 0;
			for (cutout = 0; cutout < line_distant; cutout++)
			{
				//計算斜線上的點
				if (i_all_point_map_original[show_path[i]].x >= i_all_point_map_original[show_path[j]].x && i_all_point_map_original[show_path[i]].y >= i_all_point_map_original[show_path[j]].y)
				{
					cutout_point[cutout].x = i_all_point_map_original[show_path[j]].x + (sqrt(pow((i_all_point_map_original[show_path[i]].x - i_all_point_map_original[show_path[j]].x), 2)) / (double)line_distant)*(double)cutout;
					cutout_point[cutout].y = i_all_point_map_original[show_path[j]].y + (sqrt(pow((i_all_point_map_original[show_path[i]].y - i_all_point_map_original[show_path[j]].y), 2)) / (double)line_distant)*(double)cutout;
				}

				if (i_all_point_map_original[show_path[i]].x >= i_all_point_map_original[show_path[j]].x && i_all_point_map_original[show_path[j]].y >= i_all_point_map_original[show_path[i]].y)
				{
					cutout_point[cutout].x = i_all_point_map_original[show_path[j]].x + (sqrt(pow((i_all_point_map_original[show_path[i]].x - i_all_point_map_original[show_path[j]].x), 2)) / (double)line_distant)*(double)cutout;
					cutout_point[cutout].y = i_all_point_map_original[show_path[j]].y - (sqrt(pow((i_all_point_map_original[show_path[i]].y - i_all_point_map_original[show_path[j]].y), 2)) / (double)line_distant)*(double)cutout;
				}

				if (i_all_point_map_original[show_path[j]].x >= i_all_point_map_original[show_path[i]].x && i_all_point_map_original[show_path[i]].y >= i_all_point_map_original[show_path[j]].y)
				{
					cutout_point[cutout].x = i_all_point_map_original[show_path[j]].x - (sqrt(pow((i_all_point_map_original[show_path[i]].x - i_all_point_map_original[show_path[j]].x), 2)) / (double)line_distant)*(double)cutout;
					cutout_point[cutout].y = i_all_point_map_original[show_path[j]].y + (sqrt(pow((i_all_point_map_original[show_path[i]].y - i_all_point_map_original[show_path[j]].y), 2)) / (double)line_distant)*(double)cutout;
				}

				if (i_all_point_map_original[show_path[j]].x >= i_all_point_map_original[show_path[i]].x && i_all_point_map_original[show_path[j]].y >= i_all_point_map_original[show_path[i]].y)
				{
					cutout_point[cutout].x = i_all_point_map_original[show_path[i]].x + (sqrt(pow((i_all_point_map_original[show_path[i]].x - i_all_point_map_original[show_path[j]].x), 2)) / (double)line_distant)*(double)cutout;
					cutout_point[cutout].y = i_all_point_map_original[show_path[i]].y + (sqrt(pow((i_all_point_map_original[show_path[i]].y - i_all_point_map_original[show_path[j]].y), 2)) / (double)line_distant)*(double)cutout;
				}
			}

			for (int cheak_pixel = 0; cheak_pixel < cutout; cheak_pixel++)
			{
				cheak_x = round(cutout_point[cheak_pixel].x);
				cheak_y = round(cutout_point[cheak_pixel].y);

				if ((cheak_y > i_sca_image.size() - 2) || (cheak_x > i_sca_image.size() - 2) || (cheak_y < 1 || (cheak_x < 1)))
					continue;

				if (i_sca_image[cheak_y][cheak_x] == 1 ||
					i_sca_image[cheak_y + 1][cheak_x] == 1 ||
					i_sca_image[cheak_y - 1][cheak_x] == 1 ||
					i_sca_image[cheak_y][cheak_x + 1] == 1 ||
					i_sca_image[cheak_y][cheak_x - 1] == 1 /*||
					i_sca_image[cheak_y - 1][cheak_x - 1] == 1 ||
					i_sca_image[cheak_y + 1][cheak_x + 1] == 1 ||
					i_sca_image[cheak_y - 1][cheak_x + 1] == 1 ||
					i_sca_image[cheak_y + 1][cheak_x - 1] == 1*/
					)
				{
					allpath_crash_index.push_back(1);
					crash_index = 1;
					break;
				}
				else
				{
					crash_index = 0;
				}
			}

			if (crash_index == 0)
				allpath_crash_index.push_back(0);

		}

		for (int k = allpath_crash_index.size() - 1; k >= 0; k--)  //逐一檢查路徑
		{
			if (allpath_crash_index[k] == 0)
			{
				cheak_collision_index = k + 1;  
				break;
			}
			else if (k == 0)
			{
				cheak_collision_index = 2;
			}
			else if (k == i)
			{
				cheak_collision_index = i + 1;
				break;
			}
			else
			{
				allpath_crash_index.pop_back();
			}
		}
		o_path_optimization.push_back(show_path[cheak_collision_index]);
		if (cheak_collision_index == show_path.size() - 2)
			i = cheak_collision_index;
		else
			i = cheak_collision_index;  //檢查剩餘路段是否有捷徑
		cheak_collision_index = 0;

		/*
				if (cheak_collision_index == 0) //假如都無碰撞cheak_collision_index就會等於零，直接連線終點
				{
					i = show_path.size();
					o_path_optimization.push_back(show_path[show_path.size() - 1]);
				}
				else
				{
					i = cheak_collision_index;  //檢查剩餘路段是否有捷徑
					cheak_collision_index = 0;
				}*/

	}
}


void CCollisionFreepathplanningDlg::Path_Optimization_old(vector<vector<bool>> i_sca_image, vector<CvPoint2D64f> i_all_point_map_original, int i_robot_num, vector<int>& o_path_optimization)
{


	int path_inside, break_index, line_distant;
	unsigned int cheak_x, cheak_y, cheak_collision_index = 0, final_noncollision_index = 0;
	CvPoint2D64f cutout_point[500] = { 0 };
	CvPoint2D64f temp;
	temp.x = fix_robot_end_point[i_robot_num].x / 10;
	temp.y = fix_robot_end_point[i_robot_num].y / 10;

	o_path_optimization.push_back(show_path[0]);
	show_path.push_back(i_all_point_map_original.size());
	i_all_point_map_original.push_back(temp);

	for (unsigned int i = 0; i < show_path.size(); i++)
	{
		for (unsigned int j = i + 1; j < show_path.size(); j++)
		{

			line_distant = sqrt(pow(i_all_point_map_original[show_path[i]].x - i_all_point_map_original[show_path[j]].x, 2) + pow(i_all_point_map_original[show_path[i]].y - i_all_point_map_original[show_path[j]].y, 2));

			int cutout = 0;
			for (cutout = 0; cutout < line_distant; cutout++)
			{
				//計算斜線上的點
				if (i_all_point_map_original[show_path[i]].x >= i_all_point_map_original[show_path[j]].x && i_all_point_map_original[show_path[i]].y >= i_all_point_map_original[show_path[j]].y)
				{
					cutout_point[cutout].x = i_all_point_map_original[show_path[j]].x + (sqrt(pow((i_all_point_map_original[show_path[i]].x - i_all_point_map_original[show_path[j]].x), 2)) / (double)line_distant)*(double)cutout;
					cutout_point[cutout].y = i_all_point_map_original[show_path[j]].y + (sqrt(pow((i_all_point_map_original[show_path[i]].y - i_all_point_map_original[show_path[j]].y), 2)) / (double)line_distant)*(double)cutout;
				}

				if (i_all_point_map_original[show_path[i]].x >= i_all_point_map_original[show_path[j]].x && i_all_point_map_original[show_path[j]].y >= i_all_point_map_original[show_path[i]].y)
				{
					cutout_point[cutout].x = i_all_point_map_original[show_path[j]].x + (sqrt(pow((i_all_point_map_original[show_path[i]].x - i_all_point_map_original[show_path[j]].x), 2)) / (double)line_distant)*(double)cutout;
					cutout_point[cutout].y = i_all_point_map_original[show_path[j]].y - (sqrt(pow((i_all_point_map_original[show_path[i]].y - i_all_point_map_original[show_path[j]].y), 2)) / (double)line_distant)*(double)cutout;
				}

				if (i_all_point_map_original[show_path[j]].x >= i_all_point_map_original[show_path[i]].x && i_all_point_map_original[show_path[i]].y >= i_all_point_map_original[show_path[j]].y)
				{
					cutout_point[cutout].x = i_all_point_map_original[show_path[j]].x - (sqrt(pow((i_all_point_map_original[show_path[i]].x - i_all_point_map_original[show_path[j]].x), 2)) / (double)line_distant)*(double)cutout;
					cutout_point[cutout].y = i_all_point_map_original[show_path[j]].y + (sqrt(pow((i_all_point_map_original[show_path[i]].y - i_all_point_map_original[show_path[j]].y), 2)) / (double)line_distant)*(double)cutout;
				}

				if (i_all_point_map_original[show_path[j]].x >= i_all_point_map_original[show_path[i]].x && i_all_point_map_original[show_path[j]].y >= i_all_point_map_original[show_path[i]].y)
				{
					cutout_point[cutout].x = i_all_point_map_original[show_path[i]].x + (sqrt(pow((i_all_point_map_original[show_path[i]].x - i_all_point_map_original[show_path[j]].x), 2)) / (double)line_distant)*(double)cutout;
					cutout_point[cutout].y = i_all_point_map_original[show_path[i]].y + (sqrt(pow((i_all_point_map_original[show_path[i]].y - i_all_point_map_original[show_path[j]].y), 2)) / (double)line_distant)*(double)cutout;
				}
			}

			for (int cheak_pixel = 0; cheak_pixel < cutout; cheak_pixel++)
			{
				cheak_x = round(cutout_point[cheak_pixel].x);
				cheak_y = round(cutout_point[cheak_pixel].y);

				if ((cheak_y > i_sca_image.size() - 2) || (cheak_x > i_sca_image.size() - 2) || (cheak_y < 1 || (cheak_x < 1)))
					continue;

				if (i_sca_image[cheak_y][cheak_x] == 1 ||
					i_sca_image[cheak_y + 1][cheak_x] == 1 ||
					i_sca_image[cheak_y - 1][cheak_x] == 1 ||
					i_sca_image[cheak_y][cheak_x + 1] == 1 ||
					i_sca_image[cheak_y][cheak_x - 1] == 1 ||
					i_sca_image[cheak_y - 1][cheak_x - 1] == 1 ||
					i_sca_image[cheak_y + 1][cheak_x + 1] == 1 ||
					i_sca_image[cheak_y - 1][cheak_x + 1] == 1 ||
					i_sca_image[cheak_y + 1][cheak_x - 1] == 1/**/
					)
				{
					cheak_collision_index = j - 1;  // 紀錄碰撞的上一點
					break;
				}
			}
			if (cheak_collision_index)  //紀錄此刻無碰撞捷徑後跳出搜尋
			{
				o_path_optimization.push_back(show_path[cheak_collision_index]);
				break;
			}
		}
		if (cheak_collision_index == 0) //假如都無碰撞cheak_collision_index就會等於零，直接連線終點
		{
			i = show_path.size();
			o_path_optimization.push_back(show_path[show_path.size() - 1]);
		}
		else
		{
			i = cheak_collision_index;  //檢查剩餘路段是否有捷徑
			cheak_collision_index = 0;
		}

	}
}


void CCollisionFreepathplanningDlg::MultiRobot_Path_simulation(CDC* i_pDC, IplImage * i_draw_data, vector <CPoint> i_host_path, vector<vector<bool>>  i_sca_image, vector <Point> i_save_coner, CvPoint2D64f i_robot_start_point[total_number], CvPoint2D64f i_robot_end_point[total_number], int i_car_density, vector <CPoint> &o_sim_path, vector <draw_car> &o_sim_car, IplImage *&offline_show)
{
	CvVideoWriter *master_recoder;
	char master_recoder_name[] = "master_recoder.avi";
	int FPS = 100;
	CvSize AviSize = cvSize(880, 880);
	int AviColor = 1;
	int Reentry_tag[total_number] = { 0 }, old_Reentry_tag[total_number] = { 0 };
	//	master_recoder = cvCreateVideoWriter(master_recoder_name, CV_FOURCC('D', 'I', 'V', 'X'), FPS, AviSize, AviColor);  //20180926拿掉Master

	CvVideoWriter *slave_recoder[5];
	char slave_recoder_name_1[] = "slave_recoder_1.avi";
	slave_recoder[1] = cvCreateVideoWriter(slave_recoder_name_1, CV_FOURCC('D', 'I', 'V', 'X'), FPS, AviSize, AviColor);

	char slave_recoder_name_2[] = "slave_recoder_2.avi";
	slave_recoder[2] = cvCreateVideoWriter(slave_recoder_name_2, CV_FOURCC('D', 'I', 'V', 'X'), FPS, AviSize, AviColor);

	char slave_recoder_name_3[] = "slave_recoder_3.avi";
	slave_recoder[3] = cvCreateVideoWriter(slave_recoder_name_3, CV_FOURCC('D', 'I', 'V', 'X'), FPS, AviSize, AviColor);

	char slave_recoder_name_4[] = "slave_recoder_4.avi";
	slave_recoder[4] = cvCreateVideoWriter(slave_recoder_name_4, CV_FOURCC('D', 'I', 'V', 'X'), FPS, AviSize, AviColor);

	double sampleTime = 38;
	double pixel2cm = 100;
	double rho_dot, alpha_dot, beta_dot;
	double x_here = 0, y_here = 0, zdir_here = 0, theta_here;
	double u1, u2;
	double target_pos_sim[3] = { 0 }, car_x_sim, car_y_sim, car_zdir_sim, phi_sim, rho_sim, theta_sim, alpha_sim, beta_sim;
	double vr = 0, vl = 0;
	CPoint erase_path = { 0 };
	vector <double> x_save, y_save;
	vector <draw_car>local_draw_car;
	int jump_draw = 0, state;
	//CvPoint draw_car[6];
	int carsize = 10;
	int initial_scale = 100;
	int arrive_goal_times[total_number] = { 0 };
	IplImage * live_show = NULL;

	// 	cvSetZero(draw_data);
	CvPoint draw_oringin[2];

	remove("control_pos_m_sim.txt");
	remove("car_path_efficiency.txt");
	fstream app_car_path_efficiency("car_path_efficiency.txt", ios::app);

	// 	draw_oringin[0] = cvPoint(orgin.x, orgin.y);
	// 	draw_oringin[1] = cvPoint(orgin.x, 700 - orgin.y);

// 	vector <CPoint>  jump_path_optimization_copy_sim;
// 	jump_path_optimization_copy_sim.assign(i_host_path.begin(), i_host_path.end());
// 	float x_transpos = jump_path_optimization_copy_sim[1].x;
// 	float y_transpos = jump_path_optimization_copy_sim[1].y;

// 	for (int path_num = 0; path_num < i_host_path.size() - 1; path_num++)
// 	{
	/*
			jump_draw = 0;  //20180926拿掉Master

			if (jump_path_optimization_copy_sim.size() < 2)
				break;
			else if (jump_path_optimization_copy_sim.size() < 3)
			{
				float goal_sim = atan2((jump_path_optimization_copy_sim[0].y - jump_path_optimization_copy_sim[1].y), (jump_path_optimization_copy_sim[1].x - jump_path_optimization_copy_sim[0].x));
				target_pos_sim[2] = goal_sim;  //旋轉角
			}
			else
			{
				float goal_sim = atan2((jump_path_optimization_copy_sim[1].y - jump_path_optimization_copy_sim[2].y), (jump_path_optimization_copy_sim[2].x - jump_path_optimization_copy_sim[1].x));
				target_pos_sim[2] = goal_sim;  //旋轉角

			}


			if (path_num > 0)
			{
				car_x_sim = x_here / pixel2cm * 100;
				car_y_sim = y_here / pixel2cm * 100;
				car_zdir_sim = zdir_here;
				target_pos_sim[0] = jump_path_optimization_copy_sim[1].x - jump_path_optimization_copy_sim[0].x + erase_path.x;
				target_pos_sim[1] = jump_path_optimization_copy_sim[0].y - jump_path_optimization_copy_sim[1].y + erase_path.y;

				erase_path.x = target_pos_sim[0];
				erase_path.y = target_pos_sim[1];
			}
			else
			{
				double start_pointx_sim = jump_path_optimization_copy_sim[0].x;
				double start_pointy_sim = jump_path_optimization_copy_sim[0].y;

				car_x_sim = start_pointx_sim - jump_path_optimization_copy_sim[1].x  / *+ (double)rand() / (RAND_MAX + 1.0) * 4* /;
				car_y_sim = -start_pointy_sim + jump_path_optimization_copy_sim[1].y  / *+ (double)rand() / (RAND_MAX + 1.0) * 4* /;
				car_zdir_sim = CV_PI / 2 / *+ 0.148+ (double)rand() / (RAND_MAX + 1.0) * 6* /;   //要更改初始方向角的話在這裡加   / *Camera[5] +* /
			}


			phi_sim = car_zdir_sim;
			if (phi_sim > CV_PI)		phi_sim = -2 * CV_PI + phi_sim;
			if (phi_sim < -CV_PI)		phi_sim = 2 * CV_PI + phi_sim;

			rho_sim = sqrt((car_x_sim - target_pos_sim[0])*(car_x_sim - target_pos_sim[0]) + (car_y_sim - target_pos_sim[1])*(car_y_sim - target_pos_sim[1])) *pixel2cm / 100;
			theta_sim = atan2(car_y_sim - target_pos_sim[1], car_x_sim - target_pos_sim[0]);

			alpha_sim = -phi_sim + theta_sim + CV_PI;
			if (alpha_sim > CV_PI)		alpha_sim = -2 * CV_PI + alpha_sim;
			if (alpha_sim < -CV_PI)		alpha_sim = 2 * CV_PI + alpha_sim;

			beta_sim = -(theta_sim + CV_PI) + target_pos_sim[2];
			if (beta_sim > CV_PI)		beta_sim = -2 * CV_PI + beta_sim;
			if (beta_sim < -CV_PI)		beta_sim = 2 * CV_PI + beta_sim;*/

	while (true)
	{
		if (stop_program == 1)
			break;
		/*Control_Methods(1, rho_sim, alpha_sim, beta_sim, zdir_here, vr, vl, state);  //20180926拿掉Master

		double degree_alpha = alpha_sim * 180 / CV_PI;
		double degree_beta = beta_sim * 180 / CV_PI;
		double degree_zdir_here = zdir_here * 180 / CV_PI;

		u1 = ((vr + vl) / 2);
		u2 = ((vr - vl) / 0.3);

		rho_dot = u1 * (-cos(alpha_sim));
		alpha_dot = (sin(alpha_sim) / rho_sim)*u1 - u2;
		beta_dot = -(sin(alpha_sim) / rho_sim)*u1;

		rho_sim = rho_sim + rho_dot * (float)sampleTime / 1000;
		alpha_sim = alpha_sim + alpha_dot * (float)sampleTime / 1000;
		if (alpha_sim > CV_PI)		alpha_sim = -2 * CV_PI + alpha_sim;
		if (alpha_sim < -CV_PI)		alpha_sim = 2 * CV_PI + alpha_sim;

		beta_sim = beta_sim + beta_dot * (float)sampleTime / 1000;
		if (beta_sim > CV_PI)		beta_sim = -2 * CV_PI + beta_sim;
		if (beta_sim < -CV_PI)		beta_sim = 2 * CV_PI + beta_sim;

		theta_here = -CV_PI - beta_sim + target_pos_sim[2];  //加入 target_pos_sim[2]
		if (theta_here > CV_PI)		theta_here = -2 * CV_PI + theta_here;
		if (theta_here < -CV_PI)		theta_here = 2 * CV_PI + theta_here;

		float here_scale = 1;

		x_here = cos(theta_here) * rho_sim + target_pos_sim[0];
		y_here = sin(theta_here) * rho_sim + target_pos_sim[1];
		zdir_here = -beta_sim - alpha_sim + target_pos_sim[2];  //加入 target_pos_sim[2]
		jump_draw++;
		x_save.push_back(x_here + x_transpos);
		y_save.push_back(-y_here + y_transpos);

		CPoint temp_xy;
		temp_xy.x = x_save[x_save.size() - 1];
		temp_xy.y = y_save[y_save.size() - 1];
		o_sim_path.push_back(temp_xy);

		i_robot_start_point[0].x = x_here + x_transpos;
		i_robot_start_point[0].y = -y_here + y_transpos;
		robot_zdir[0] = zdir_here;*/

		for (int robot_number = 1; robot_number < total_number; robot_number++)
		{

			if (abs(i_robot_start_point[robot_number].x - i_robot_end_point[robot_number].x) > 2 || abs(i_robot_start_point[robot_number].y - i_robot_end_point[robot_number].y) > 2)
			{
				servant_path(robot_number, i_draw_data, i_sca_image, i_save_coner, i_robot_start_point, i_robot_end_point, o_sim_path, carsize, slave_recoder[robot_number], arrive_goal_times);
			}
			else  //可以讓機器來回行走，交換起點與終點
			{
				Reentry_tag[robot_number] = all_robot_path_double[robot_number].size();  //紀錄折返的當下資料位置
				car_path_efficiency(fix_robot_start_point[robot_number], i_robot_end_point[robot_number], Reentry_tag[robot_number], old_Reentry_tag[robot_number], robot_number);

				old_Reentry_tag[robot_number] = Reentry_tag[robot_number];

				arrive_goal_times[robot_number]++;
				memcpy(&i_robot_end_point[robot_number], &fix_robot_start_point[robot_number], sizeof(fix_robot_start_point[robot_number]));
				memcpy(&fix_robot_start_point[robot_number], &fix_robot_end_point[robot_number], sizeof(fix_robot_end_point[robot_number]));
				memcpy(&fix_robot_end_point[robot_number], &i_robot_end_point[robot_number], sizeof(i_robot_end_point[robot_number]));
				servant_path(robot_number, i_draw_data, i_sca_image, i_save_coner, i_robot_start_point, i_robot_end_point, o_sim_path, carsize, slave_recoder[robot_number], arrive_goal_times);

				app_car_path_efficiency << robot_number << " " << arrive_goal_times[robot_number] << " " << show_path_efficiency[robot_number] << endl;

			}


		}

		if (arrive_goal_times[1] > 50)
			break;



		// 			app_pos_m_sim
		// 				<< zdir_here << "  "
		// 				<< x_here << "  "
		// 				<< y_here << "  "
		// 				<< endl;

/*
		draw_car draw_the_car[total_number];

		//if (jump_draw % i_car_density == 1)  //i_car_density 顯示頻率，不顯示就設很大
		if (0)
		{
			live_show = cvCloneImage(offline_show);

			for (int k = 0; k < total_number; k++)
				simulation_car(live_show, i_robot_start_point, robot_zdir, k, carsize, arrive_goal_times);

			local_draw_car.push_back(draw_the_car[0]);

			CvvImage show1;
			show1.CopyOf(live_show);
			show1.Show(*i_pDC, 0, 0, live_show->width, live_show->height);

			//	cvWriteFrame(master_recoder, live_show);

			cvReleaseImageData(live_show);
		}*/
		// 			if ((abs(x_here - target_pos_sim[0] * pixel2cm / 100) < 1 && abs(y_here - target_pos_sim[1] * pixel2cm / 100) < 1) /*|| jump_draw > 50 || zdir_here < 0.2*/)  //20180926拿掉Master
		// 			{
		// 				jump_path_optimization_copy_sim.erase(jump_path_optimization_copy_sim.begin(), jump_path_optimization_copy_sim.begin() + 1);
		// 				break;
		// 			}
	}
	//static_draw_car.assign(local_draw_car.begin(), local_draw_car.end());

// 
// 		for (int i = 0; i < x_save.size() - 1; i++)
// 		{
// 			cvLine(offline_show, cvPoint((int)x_save[i], (int)y_save[i]), cvPoint((int)x_save[i + 1], (int)y_save[i + 1]), CV_RGB(150, 150, 255), 2);
// 		}

//	}

	app_car_path_efficiency.close();

	//	cvReleaseVideoWriter(&master_recoder);
	cvReleaseVideoWriter(&slave_recoder[0]);
	cvReleaseVideoWriter(&slave_recoder[1]);
	cvReleaseVideoWriter(&slave_recoder[2]);
	cvReleaseVideoWriter(&slave_recoder[3]);
	cvReleaseVideoWriter(&slave_recoder[4]);
	cvReleaseImage(&live_show);
	o_sim_car = local_draw_car;
	local_draw_car.clear();
	return;

}

void CCollisionFreepathplanningDlg::ServantRobot_Path_simulation(vector <CPoint> i_Servant_path, int i_robot_num, CvPoint2D64f &o_ServantRobot_pos, CvPoint2D64f i_robot_start_point, double &io_zdir)
{
	//	remove("子機器人輸出.txt");
//	fstream app_ServantRobot("子機器人輸出.txt", ios::app);

	double sampleTime = 40;
	double pixel2cm = 100;
	double rho_dot, alpha_dot, beta_dot;
	double x_here, y_here, zdir_here = 0, theta_here;
	double u1, u2;
	double target_pos_sim[3] = { 0 }, car_x_sim, car_y_sim, car_zdir_sim, phi_sim, rho_sim, theta_sim, alpha_sim, beta_sim;
	double vr = 0, vl = 0;
	double x_save, y_save;
	vector <draw_car>local_draw_car;
	int jump_draw = 0, state;
	//CvPoint draw_car[6];
	int carsize = 20;
	int initial_scale = 100;
	IplImage * live_show = NULL;

	// 	cvSetZero(draw_data);

	// 	remove("control_pos_m_sim.txt");
	// 	fstream app_pos_m_sim("control_pos_m_sim.txt", ios::app);

	// 	draw_oringin[0] = cvPoint(orgin.x, orgin.y);
	// 	draw_oringin[1] = cvPoint(orgin.x, 700 - orgin.y);

	vector <CPoint>  jump_path_optimization_copy_sim;
	jump_path_optimization_copy_sim.assign(i_Servant_path.begin(), i_Servant_path.end());


	double x_transpos = jump_path_optimization_copy_sim[1].x;
	double y_transpos = jump_path_optimization_copy_sim[1].y;


	double distant_check = sqrt((i_robot_start_point.x - jump_path_optimization_copy_sim[1].x - target_pos_sim[0])*(i_robot_start_point.x - jump_path_optimization_copy_sim[1].x - target_pos_sim[0]) + (-i_robot_start_point.y + jump_path_optimization_copy_sim[1].y - target_pos_sim[1])*(-i_robot_start_point.y + jump_path_optimization_copy_sim[1].y - target_pos_sim[1]));



	if (jump_path_optimization_copy_sim.size() < 3 || distant_check > 50)
	{
		double goal_sim = atan2((jump_path_optimization_copy_sim[0].y - jump_path_optimization_copy_sim[1].y), (jump_path_optimization_copy_sim[1].x - jump_path_optimization_copy_sim[0].x));
		target_pos_sim[2] = goal_sim;  //旋轉角
	}
	else
	{
		double goal_sim = atan2((jump_path_optimization_copy_sim[1].y - jump_path_optimization_copy_sim[2].y), (jump_path_optimization_copy_sim[2].x - jump_path_optimization_copy_sim[1].x));
		target_pos_sim[2] = goal_sim;  //旋轉角

	}


	double start_pointx_sim = i_robot_start_point.x;
	double start_pointy_sim = i_robot_start_point.y;

	car_x_sim = start_pointx_sim - jump_path_optimization_copy_sim[1].x  /*+ (double)rand() / (RAND_MAX + 1.0) * 4*/;
	car_y_sim = -start_pointy_sim + jump_path_optimization_copy_sim[1].y  /*+ (double)rand() / (RAND_MAX + 1.0) * 4*/;
	car_zdir_sim = io_zdir; /*+ 0.148+ (double)rand() / (RAND_MAX + 1.0) * 6*/;   //要更改初始方向角的話在這裡加   /*Camera[5] +*/

	if (distant_check > 50)
	{
		target_pos_sim[0] = 5 * car_x_sim / 6;
		target_pos_sim[1] = 5 * car_y_sim / 6;
	}

	phi_sim = car_zdir_sim;
	if (phi_sim > CV_PI)		phi_sim = -2 * CV_PI + phi_sim;
	if (phi_sim < -CV_PI)		phi_sim = 2 * CV_PI + phi_sim;

	rho_sim = sqrt((car_x_sim - target_pos_sim[0])*(car_x_sim - target_pos_sim[0]) + (car_y_sim - target_pos_sim[1])*(car_y_sim - target_pos_sim[1])) *pixel2cm / 100;
	theta_sim = atan2(car_y_sim - target_pos_sim[1], car_x_sim - target_pos_sim[0]);

	alpha_sim = -phi_sim + theta_sim + CV_PI;
	if (alpha_sim > CV_PI)		alpha_sim = -2 * CV_PI + alpha_sim;
	if (alpha_sim < -CV_PI)		alpha_sim = 2 * CV_PI + alpha_sim;

	beta_sim = -(theta_sim + CV_PI) + target_pos_sim[2];
	if (beta_sim > CV_PI)		beta_sim = -2 * CV_PI + beta_sim;
	if (beta_sim < -CV_PI)		beta_sim = 2 * CV_PI + beta_sim;

	//	app_ServantRobot << rho_sim << " " << alpha_sim << " " << beta_sim << " " << car_x_sim << " " << car_y_sim << " " << car_zdir_sim;
		//----------------------------------------------------------------

	Control_Methods(2, rho_sim, alpha_sim, beta_sim, zdir_here, vr, vl, state);

	for (int j = i_robot_num; j < total_number; j++) //如果高優先權機器人前面有較低優先權之機器人，則放慢速度或直接暫停前進
	{
		if (car_front_num[i_robot_num][j] == 1)
		{
			vr = 0;
			vl = 0;
		}
	}



	double degree_alpha = alpha_sim * 180 / CV_PI;
	double degree_beta = beta_sim * 180 / CV_PI;
	double degree_zdir_here = zdir_here * 180 / CV_PI;

	u1 = ((vr + vl) / 2);
	u2 = ((vr - vl) / 0.3);

	rho_dot = u1 * (-cos(alpha_sim));
	alpha_dot = (sin(alpha_sim) / rho_sim)*u1 - u2;
	beta_dot = -(sin(alpha_sim) / rho_sim)*u1;

	rho_sim = rho_sim + rho_dot * sampleTime / 1000.0;
	alpha_sim = alpha_sim + alpha_dot * sampleTime / 1000.0;
	if (alpha_sim > CV_PI)		alpha_sim = -2 * CV_PI + alpha_sim;
	if (alpha_sim < -CV_PI)		alpha_sim = 2 * CV_PI + alpha_sim;

	beta_sim = beta_sim + beta_dot * sampleTime / 1000.0;
	if (beta_sim > CV_PI)		beta_sim = -2 * CV_PI + beta_sim;
	if (beta_sim < -CV_PI)		beta_sim = 2 * CV_PI + beta_sim;

	theta_here = -CV_PI - beta_sim + target_pos_sim[2];  //加入 target_pos_sim[2]
	if (theta_here > CV_PI)		theta_here = -2 * CV_PI + theta_here;
	if (theta_here < -CV_PI)		theta_here = 2 * CV_PI + theta_here;

	float here_scale = 1;

	x_here = cos(theta_here) * rho_sim + target_pos_sim[0];
	y_here = sin(theta_here) * rho_sim + target_pos_sim[1];
	zdir_here = -beta_sim - alpha_sim + target_pos_sim[2];  //加入 target_pos_sim[2]
	io_zdir = zdir_here;

	//	app_ServantRobot << "    " << rho_sim << " " << alpha_sim << " " << beta_sim << " " << x_here << " " << y_here << " " << zdir_here << endl;

	o_ServantRobot_pos.x = (x_here + x_transpos);
	o_ServantRobot_pos.y = (-y_here + y_transpos);

	//	app_ServantRobot.close();
}

void CCollisionFreepathplanningDlg::servant_path(int i_robot_num, IplImage * i_pGrayImg, vector<vector<bool>> i_sca_image, vector <Point> i_save_coner, CvPoint2D64f i_robot_start_point[total_number], CvPoint2D64f i_robot_end_point[total_number], vector <CPoint> i_sim_path, int i_carsize, CvVideoWriter *i_slave_recoder, int i_arrive_goal_times[total_number])
{
	int corner_count = 0;
	int line_count = 0;  //有多少VD線段
	int new_input_index = 0;
	int path_optimization_size_change = 0;
	char num1[50];
	double Data[8000] = { 0 };
	CElement* pElement = 0;
	CvPoint2D32f *point2 = 0;
	CvPoint2D64f savepoint1[3000] = { 0.0 }, savepoint2[3000] = { 0.0 };
	CvPoint2D64f new_savepoint1[3000] = { 0.0 }, new_savepoint2[3000] = { 0.0 };
	vector <CPoint> all_point_map;
	vector <CPoint>  jump_path_optimization, jump_path_optimization_simulation;
	vector <Point> save_coner;
	vector <CvPoint2D64f> all_point_map_original;
	vector <int> path_optimization;
	IplImage * draw_data6 = NULL;
	IplImage * draw_data5 = NULL;
	IplImage * draw_data4 = NULL;
	IplImage * draw_data3 = NULL;
	IplImage * draw_data2 = NULL;
	IplImage * draw_data1 = NULL;
	CvFont Font1 = cvFont(2, 1);

	//	draw_data3 = cvCreateImage(cvGetSize(i_pGrayImg), 8, 3);
	draw_data2 = cvCreateImage(cvSize(880, 880), 8, 3);
	draw_data4 = cvCreateImage(cvSize(88, 88), 8, 3);
	draw_data5 = cvCreateImage(cvSize(440, 440), 8, 3);
	draw_data6 = cvCreateImage(cvSize(88, 88), 8, 3);
	draw_data1 = cvCreateImage(cvSize(88, 88), 8, 1);
	draw_data3 = cvCloneImage(i_pGrayImg);

	Point host_position[50];
	Point other_position;

	//----------------------------清除數據----------------------------------------------------------------------
	memset(d, 0, sizeof(d));
	memset(parent, 0, sizeof(parent));
	memset(visit, 0, sizeof(visit));
	memset(w, 0, sizeof(w));
	show_path.clear();
	point_data_.clear();
	vd.clear();
	//-------------------------------------------------------------------------------------------------------------




/*
	host_position[0].x = i_sim_path[i_sim_path.size() - 1].x / 10;  //提取Master機器人現在位置
	host_position[0].y = i_sim_path[i_sim_path.size() - 1].y / 10;  //20180926拿掉Master

		for (int i = 0; i < 2; i++)
		{
			for (int j = -2; j < 2; j++)
			{
				host_position[i * 4 + (j + 3)].x = host_position[0].x + j;
				host_position[i * 4 + (j + 3)].y = host_position[0].y - i;
			}
		}

		//插入主機器人之位置與障礙
		for (int i = 0; i < 8; i++)
		{
			//		save_coner.push_back(host_position[i]);
			i_sca_image[host_position[i].y][host_position[i].x] = 1;
	 */

	delay_own_obstacle[i_robot_num][1] = i_robot_start_point[i_robot_num];
	double own_dis = sqrt(pow(delay_own_obstacle[i_robot_num][1].x - delay_own_obstacle[i_robot_num][0].x, 2) + pow(delay_own_obstacle[i_robot_num][1].y - delay_own_obstacle[i_robot_num][0].y, 2));

	//新增延遲更新系統，當自己的coner離自己一段距離後才更新到當下位置
	//綠明顯---------------------------------------------
	//綠明顯~~~~~~~~~~~~~~~~~~~~~~~~~~
	//綠明顯~~~~~~~~~~~~~~~~~~~~~~~~~~
	if (own_dis > 30)
		delay_own_obstacle[i_robot_num][0] = i_robot_start_point[i_robot_num];
	//綠明顯~~~~~~~~~~~~~~~~~~~~~~~~~~
	//綠明顯~~~~~~~~~~~~~~~~~~~~~~~~~~
	//綠明顯~~~~~~~~~~~~~~~~~~~~~~~~~~
	//綠明顯~~~~~~~~~~~~~~~~~~~~~~~~~~

	for (int i = i_robot_num; i > 0; i--)  //新增主從機器人之障礙格點
	{
		if (delay_own_obstacle[i][0].x == 0 && delay_own_obstacle[i][0].y == 0 || i == i_robot_num)
			continue;

		other_position.x = delay_own_obstacle[i][0].x / 10;
		other_position.y = delay_own_obstacle[i][0].y / 10;

		save_coner.push_back(other_position);
		i_sca_image[other_position.y][other_position.x] = 1;
		i_sca_image[other_position.y + 1][other_position.x + 1] = 1;
		i_sca_image[other_position.y + 1][other_position.x] = 1;
		i_sca_image[other_position.y][other_position.x + 1] = 1;
		i_sca_image[other_position.y - 1][other_position.x - 1] = 1;
		i_sca_image[other_position.y - 1][other_position.x] = 1;
		i_sca_image[other_position.y][other_position.x - 1] = 1;
		i_sca_image[other_position.y - 1][other_position.x + 1] = 1;
		i_sca_image[other_position.y + 1][other_position.x - 1] = 1;
		cvSetReal2D(draw_data3, other_position.y, other_position.x, 0);
		cvSetReal2D(draw_data3, other_position.y + 1, other_position.x + 1, 0);
		cvSetReal2D(draw_data3, other_position.y + 1, other_position.x, 0);
		cvSetReal2D(draw_data3, other_position.y, other_position.x + 1, 0);
		cvSetReal2D(draw_data3, other_position.y - 1, other_position.x - 1, 0);
		cvSetReal2D(draw_data3, other_position.y - 1, other_position.x, 0);
		cvSetReal2D(draw_data3, other_position.y, other_position.x - 1, 0);
		cvSetReal2D(draw_data3, other_position.y + 1, other_position.x - 1, 0);
		cvSetReal2D(draw_data3, other_position.y - 1, other_position.x + 1, 0);
	}

	// 	other_position.x = i_robot_start_point[i_robot_num].x / 10;
	// 	other_position.y = i_robot_start_point[i_robot_num].y / 10;
	// 	save_coner.push_back(other_position);

		//輸入二值資料，輸出角點
	find_coner(i_sca_image, save_coner, 4);
	//插入其他子機器人之位置與障礙
	// 	save_coner.push_back(host_position[0]);     //新增Master機器人位置當障礙物
	// 	i_sca_image[host_position[0].y][host_position[0].x] = 1;  //20180926拿掉Master
	// 	cvSetReal2D(draw_data3, host_position[0].y, host_position[0].x, 0);

	//------------------------------------------------------
	other_position.x = delay_own_obstacle[i_robot_num][0].x / 10; //加入自身當障礙物創造路徑
	other_position.y = delay_own_obstacle[i_robot_num][0].y / 10;
	save_coner.push_back(other_position);
	i_sca_image[other_position.y][other_position.x] = 1;
	cvSetReal2D(draw_data3, other_position.y, other_position.x, 0);
	//------------------------------------------------------

	//將角點轉換為準備要丟入Voronoi運算的格式
	trans2Voronoi(i_sca_image, save_coner, Data, 8);
	//計算狹義Voronoi，輸入角點資料與邊界，輸出兩個矩陣
//	Voronoi_calculate(Data, 88, 88, savepoint1, savepoint2, line_count);
	//計算狹義boost_Voronoi，輸入角點資料與邊界，輸出兩個矩陣，新版使用boost c++
	boost_Voronoi_calculate(88, 88, savepoint1, savepoint2, line_count);
	//計算廣義Voronoi，待改
	Generalized_Voronoi(i_sca_image, savepoint1, savepoint2, line_count, new_input_index, new_savepoint1, new_savepoint2);
	//VD點會破碎，將其重新聚合
	Match_point(line_count, new_input_index, new_savepoint1, new_savepoint2, 1);
	//Dijkstra路徑搜尋，輸入點連接資訊跟數量
	Dijkstra_path_planning(i_robot_num, i_robot_start_point, i_robot_end_point, new_savepoint1, new_savepoint2, new_input_index, all_point_map, all_point_map_original);

	//------------------------------------------------------
	save_coner.pop_back();  //再刪除自身障礙物來連線
	i_sca_image[other_position.y][other_position.x] = 0;
	cvSetReal2D(draw_data3, other_position.y, other_position.x, 1);
	//------------------------------------------------------

	//路徑優化，輸入二值資訊與原本路徑
	//綠明顯~~~~~~~~~~~~~~~~~~~~~~~~~~
	//綠明顯~~~~~~~~~~~~~~~~~~~~~~~~~~
	//綠明顯~~~~~~~~~~~~~~~~~~~~~~~~~~
	Path_Optimization_20181127(i_sca_image, all_point_map_original, i_robot_num, path_optimization); //新版捷徑搜尋，可以跳過障礙物之後繼續搜尋路徑
//	Path_Optimization_old(i_sca_image, all_point_map_original, i_robot_num, path_optimization);
	//綠明顯~~~~~~~~~~~~~~~~~~~~~~~~~~
	//綠明顯~~~~~~~~~~~~~~~~~~~~~~~~~~
	//綠明顯~~~~~~~~~~~~~~~~~~~~~~~~~~

	IplConvKernel *pKernel_small = NULL;
	//	pKernel_small = cvCreateStructuringElementEx(11, 11, 6, 6, CV_SHAPE_RECT, NULL);
	pKernel_small = cvCreateStructuringElementEx(3, 3, 1, 1, CV_SHAPE_RECT, NULL);

	// 	cvResize(draw_data3, draw_data1, CV_INTER_NN);
	// 	cvCvtColor(draw_data1, draw_data4, CV_GRAY2RGB);
	// 	cvDilate(draw_data4, draw_data2, pKernel_small, 2);
	// 	cvAddWeighted(draw_data2, 0.5, draw_data4, 0.5, 0, draw_data2);


	cvDilate(draw_data3, draw_data1, pKernel_small, 2);
	cvCvtColor(draw_data1, draw_data4, CV_GRAY2RGB);
	cvCvtColor(draw_data3, draw_data6, CV_GRAY2RGB);
	cvAddWeighted(draw_data4, 0.5, draw_data6, 0.5, 0, draw_data4);
	cvResize(draw_data4, draw_data2, CV_INTER_NN);
	//-------------------------------------------繪圖---------------------------------------
	CPoint final_goal_point;
	final_goal_point.x = fix_robot_end_point[i_robot_num].x;
	final_goal_point.y = fix_robot_end_point[i_robot_num].y;
	all_point_map.push_back(final_goal_point);
	// 	for (int i = 0; i < save_coner.size(); i++)  //角點圖
	// 	{
	// 		cvLine(draw_data2, cvPoint(save_coner[i].x * 10, save_coner[i].y * 10), cvPoint(save_coner[i].x * 10, save_coner[i].y * 10), CV_RGB(0, 250, 250), 8);
	// 	}
	for (int i = 0; i < line_count; i++)   //VD圖
	{
		cvLine(draw_data2, cvPoint(savepoint1[i].x * 10, savepoint1[i].y * 10), cvPoint(savepoint2[i].x * 10, savepoint2[i].y * 10), CV_RGB(0, 0, 255), 1);
	}
	for (int i = 0; i < new_input_index; i++)  //GVD圖
	{
		cvLine(draw_data2, cvPoint(new_savepoint1[i].x * 10, new_savepoint1[i].y * 10), cvPoint(new_savepoint2[i].x * 10, new_savepoint2[i].y * 10), CV_RGB(250, 100, 100), 1);
	}

	for (unsigned int path_index = 0; path_index < show_path.size() - 1; path_index++) //畫出路徑圖
	{
		cvLine(draw_data2, cvPoint(all_point_map[show_path[path_index]].x, all_point_map[show_path[path_index]].y), cvPoint(all_point_map[show_path[path_index + 1]].x, all_point_map[show_path[path_index + 1]].y), CV_RGB(0, 0, 255), 3);
	}




	for (unsigned int path_opt = 0; path_opt < path_optimization.size(); path_opt++) //路徑優化顯示
	{
		jump_path_optimization.push_back(all_point_map[path_optimization[path_opt]]);

		if (path_opt == path_optimization.size() - 1)
			break;
		else
		{
			cvLine(draw_data2, cvPoint(all_point_map[path_optimization[path_opt]].x, all_point_map[path_optimization[path_opt]].y), cvPoint(all_point_map[path_optimization[path_opt + 1]].x, all_point_map[path_optimization[path_opt + 1]].y), CV_RGB(0, 150, 0), 3);
		}
	}
	/*
		float pathpoint_dis = 200;
		int jump_num2 = 20000, jump_num;
		int path_opt = 0;

		for (int path_opt = 0; path_opt < path_optimization.size(); path_opt++) //路徑優化搜尋，20181004不知幹嘛 故拔掉
		{
			if (path_optimization_size_change != path_optimization.size())
			{
				jump_num = 3000000;
			}

			if (pathpoint_dis > 5)  //16
			{
				//					first_in = sqrtf(pow((all_point_map[path_optimization[0]].x - all_point_map[path_optimization[1]].x), 2) + pow((all_point_map[path_optimization[0]].y - all_point_map[path_optimization[1]].y), 2)) / 2;
				jump_num2 = path_optimization[path_opt];

				if (jump_num2 != jump_num)
					jump_path_optimization.push_back(all_point_map[path_optimization[path_opt]]);

			}
			else
			{
				jump_num = path_optimization[1];
				//					jump_path_optimization_save = all_point_map[path_optimization[1]];

			}
			if (path_opt == path_optimization.size() - 1)
				break;
			else
			{
				pathpoint_dis = sqrtf(pow((all_point_map[path_optimization[0]].x - all_point_map[path_optimization[path_opt + 1]].x), 2) + pow((all_point_map[path_optimization[0]].y - all_point_map[path_optimization[path_opt + 1]].y), 2));


				cvLine(draw_data2, cvPoint(all_point_map[path_optimization[path_opt]].x, all_point_map[path_optimization[path_opt]].y), cvPoint(all_point_map[path_optimization[path_opt + 1]].x, all_point_map[path_optimization[path_opt + 1]].y), CV_RGB(0, 150, 0), 3);
				//				cvCircle(draw_data2, cvPoint(all_point_map[path_optimization[path_opt]].x, all_point_map[path_optimization[path_opt]].y), 7, CV_RGB(0, 150, 0), 2);

			}
		}
		*/
	if (jump_path_optimization.size() < 2)
	{
		jump_path_optimization.push_back(CPoint(i_robot_end_point[i_robot_num].x, i_robot_end_point[i_robot_num].y));
		// 		cvReleaseImage(&draw_data2);  //原為到達點之後就停止，已改為無限循環
		// 		cvReleaseImage(&draw_data1);
		// 		cvReleaseImage(&draw_data3);
		// 		return;
	}

	CvPoint2D64f temp_xy_double;
	CvPoint temp_xy;
	ServantRobot_Path_simulation(jump_path_optimization, i_robot_num, temp_xy_double, i_robot_start_point[i_robot_num], robot_zdir[i_robot_num]);

	i_robot_start_point[i_robot_num] = temp_xy_double;
	temp_xy.x = i_robot_start_point[i_robot_num].x;
	temp_xy.y = i_robot_start_point[i_robot_num].y;
	all_robot_path[i_robot_num].push_back(temp_xy);
	all_robot_dir[i_robot_num].push_back(robot_zdir[i_robot_num]);
	all_robot_path_double[i_robot_num].push_back(temp_xy_double);

	for (int i = 0; i < total_number; i++)
		simulation_car(draw_data2, i_robot_start_point, robot_zdir, i, i_carsize, i_arrive_goal_times);

	car_crash_detector(i_robot_start_point, i_carsize, car_crash_times, crash_num);
	sprintf_s(num1, "Crash times %d", car_crash_times);
	cvPutText(draw_data2, num1, cvPoint(30, 30), &Font1, Scalar(50, 50, 255));

	if (i_robot_num == 2)
	{
		cvShowImage("draw_data", draw_data2); // 動態顯示模擬
		cvWaitKey(5); // 停留視窗
	}

	//	cvResize(draw_data2, draw_data5, CV_INTER_LINEAR);
	//	cvWriteFrame(i_slave_recoder, draw_data2);  //存成影片

	cvReleaseImage(&draw_data1);
	cvReleaseImage(&draw_data2);
	cvReleaseImage(&draw_data3);
	cvReleaseImage(&draw_data4);
	cvReleaseImage(&draw_data5);
	cvReleaseImage(&draw_data6);
}

void CCollisionFreepathplanningDlg::Control_Methods(int control_type, double i_rho, double i_alpha, double i_beta, double i_phi, double &o_vr, double &o_vl, int &o_state)
{
	//----------------------第二模式新增之判斷------------------------------------------------------------
	float P2_1 = 1.1737;
	float P2_2 = 1.4317;
	float P2_3 = 0.2422;
	float P2_4 = 0.4223;

	float Vt2 = i_rho * i_rho * P2_1 +
		i_alpha * i_alpha*P2_2 +
		i_alpha * i_phi*P2_4 +
		i_alpha * i_phi*P2_4 +
		i_phi * i_phi*P2_3;

	int rho_gain;
	int times = 0;
	float M1, M2, N1, N2;
	float alpha_gain;
	float beta_gain;
	//-----------------------------------------------------------------------------------------------------


	if (control_type == 1)
	{
		//原始線性控制
		o_vr = 2 * i_rho + 0.15 * (20 * i_alpha - 2 * (i_beta));
		o_vl = 2 * i_rho - 0.15 * (20 * i_alpha - 2 * (i_beta));
		//		o_vr = o_vr / 4;
		//		o_vl = o_vl / 4;


		//-------------------------normalize----------------------------------
		double temp = abs(o_vr) + abs(o_vl);
		o_vr = 50 * o_vr / temp;
		o_vl = 50 * o_vl / temp;
	}
	else if (control_type == 2)  //超簡單控制模式
	{
		if (abs(i_alpha) > 0.03)
		{
			o_vr = 0 * i_rho + 0.15 * (1 * i_alpha - 0 * (i_beta));
			o_vl = 0 * i_rho - 0.15 * (1 * i_alpha - 0 * (i_beta));
			double temp = abs(o_vr) + abs(o_vl);
			o_vr = 0.4 * o_vr / temp;
			o_vl = 0.4 * o_vl / temp;
		}
		else
		{
			o_vr = 1 * i_rho + 0.15 * (0.01 * i_alpha - 0 * (i_beta));
			o_vl = 1 * i_rho - 0.15 * (0.01 * i_alpha - 0 * (i_beta));
			double temp = abs(o_vr) + abs(o_vl);
			o_vr = 150 * o_vr / temp;
			o_vl = 150 * o_vl / temp;
		}
	}
	else if (control_type == 3)
	{
		//判斷要切換哪種模式TS-FUZZY
		if (times == 0 || i_alpha > (CV_PI / 10) || i_alpha < (-CV_PI / 10))
		{
			if (o_state == 3)
				o_state = 3;

			else if (o_state == 2)
				o_state = 2;

			else
				o_state = 1;
		}
		else if ((i_rho > 1 || Vt2 > 4))  //gamma = 2
		{
			if (o_state == 3)
				o_state = 3;
			else
				o_state = 2;
		}
		else
			o_state = 3;

		if (i_alpha < CV_PI / 10 && i_alpha>-CV_PI / 10)
			times = 1; //一開始必定Mode1，其餘只看角度來決定是否進入Mode1

					   //切換式TS-Fuzzy控制
		alpha_gain = 1.9161;

		switch (o_state)
		{
		case 1:  //旋轉Mode
				 //		alpha_gain = 10 * (i_alpha / pi);

			o_vr = 0 * i_rho + 0.15 * (alpha_gain*i_alpha - 0 * (i_beta));
			o_vl = 0 * i_rho - 0.15 * (alpha_gain*i_alpha - 0 * (i_beta));
			o_vr = o_vr * 1.4;
			o_vl = o_vl * 1.4;

			break;

		case 2:  //直線Mode
				 // 			rho_gain = (4 * i_rho / 200) + (1 - (i_rho / 200));
				 // 			if (i_rho > 150) rho_gain = 4;
			rho_gain = 1.0331;
			o_vr = rho_gain * i_rho + 0.15 * (alpha_gain*i_alpha - 0 * (i_beta));
			o_vl = rho_gain * i_rho - 0.15 * (alpha_gain*i_alpha - 0 * (i_beta));
			o_vr = o_vr * 0.1;
			o_vl = o_vl * 0.1;
			break;

		case 3:  //PDC Mode
				 //			rho_gain = (3 * i_rho / 125) + (1 - (i_rho / 125));
			rho_gain = 1.0331;
			M1 = (cos(i_alpha) - 0.031415926) / (1 - 0.031415926);
			M2 = 1 - M1;
			//		M2 = (1 - cos(i_alpha)) / (1 - 0.031415926);

			if (i_alpha == 0)
				N1 = 1;
			else
				N1 = (0.49*CV_PI * sin(i_alpha) - sin(0.49*CV_PI)*i_alpha) / (i_alpha*(0.49*CV_PI - sin(0.49*CV_PI)));

			N2 = 1 - N1;

			// 		alpha_gain = M1*N1*5.8766 +
			// 			                  M1*N2*5.6385 +
			// 			                  M2*N1*5.8766 +
			// 			                  M2*N2*5.6385;
			// 
			// 		beta_gain = M1*N1*1.1052 +
			// 			                M1*N2*1.0776 +
			// 			                M2*N1*1.1052 +
			// 			                M2*N2*1.0776;

			alpha_gain = M1 * N1*1.2833 +
				M1 * N2*1.1022 +
				M2 * N1* 1.2833 +
				M2 * N2*1.1022;

			beta_gain = -M1 * N1*0.0487 +
				-M1 * N2*0.0517 +
				-M2 * N1*0.0487 +
				-M2 * N2*0.0517;

			// 			beta_gain = -M1*N1*1.2833 +
			// 				-M1*N2*1.1022 +
			// 				-M2*N1* 1.2833 +
			// 				-M2*N2*1.1022;
			// 
			// 			alpha_gain = M1*N1*0.0487 +
			// 				M1*N2*0.0517 +
			// 				M2*N1*0.0487 +
			// 				M2*N2*0.0517;


			o_vr = rho_gain * i_rho + 0.15 * (alpha_gain * i_alpha - beta_gain * (i_beta));
			o_vl = rho_gain * i_rho - 0.15 * (alpha_gain * i_alpha - beta_gain * (i_beta));
			o_vr = o_vr * 0.3;
			o_vl = o_vl * 0.3;
			break;

		default:
			break;
		}
	}

}

void CCollisionFreepathplanningDlg::simulation_car(IplImage *& live_show, CvPoint2D64f i_robot_start_point[total_number], double i_robot_zdir[total_number], int i_car_num, int i_carsize, int i_arrive_goal_times[total_number])
{
	//	fstream app_pos_m_sim("control_pos_m_sim.txt", ios::app);
	int total_arrive_times = 0;
	CvPoint car[6];
	char num1[200];
	CvFont Font1 = cvFont(2, 2);
	CvFont Font2 = cvFont(1, 1);
	CvPoint temp_path, detect_area[total_number][2] = { 0 };
	car_front_detector(i_robot_start_point, robot_zdir, 1, i_carsize * 6, detect_area);  //偵測前方是誰，並輸出繪圖用偵測範圍

	car[0] = cvPoint(i_robot_start_point[i_car_num].x + i_carsize * cos(robot_zdir[i_car_num] + 0.7854), i_robot_start_point[i_car_num].y - i_carsize * sin(robot_zdir[i_car_num] + 0.7854));
	car[1] = cvPoint(i_robot_start_point[i_car_num].x + i_carsize * cos(0.7854 - robot_zdir[i_car_num]), i_robot_start_point[i_car_num].y + i_carsize * sin(0.7854 - robot_zdir[i_car_num]));
	car[2] = cvPoint(i_robot_start_point[i_car_num].x - i_carsize * cos(robot_zdir[i_car_num] + 0.7854), i_robot_start_point[i_car_num].y + i_carsize * sin(robot_zdir[i_car_num] + 0.7854));
	car[3] = cvPoint(i_robot_start_point[i_car_num].x - i_carsize * cos(0.7854 - robot_zdir[i_car_num]), i_robot_start_point[i_car_num].y - i_carsize * sin(0.7854 - robot_zdir[i_car_num]));
	car[4] = cvPoint(i_robot_start_point[i_car_num].x, i_robot_start_point[i_car_num].y);
	car[5] = cvPoint(i_robot_start_point[i_car_num].x + 20 * cos(robot_zdir[i_car_num]), i_robot_start_point[i_car_num].y - 20 * sin(robot_zdir[i_car_num]));


	/*if (i_car_num != 3)
		app_pos_m_sim << car[4].x << " " << car[4].y << " " << robot_zdir[i_car_num] << " ";
	else
		app_pos_m_sim << car[4].x << " " << car[4].y << " " << robot_zdir[i_car_num] << endl;*/

		// 	temp_path.x = i_robot_start_point[i_car_num].x;
		// 	temp_path.y = i_robot_start_point[i_car_num].y;
		// 
		// 	all_robot_path[i_car_num].push_back(temp_path);
		// 	all_robot_dir[i_car_num].push_back(robot_zdir[i_car_num]);

	if (i_car_num == 0)
	{
		return;
		/*
				cvLine(live_show, car[0], car[1], CV_RGB(0, 255, 0), 2);
				cvLine(live_show, car[1], car[2], CV_RGB(0, 255, 0), 2);
				cvLine(live_show, car[2], car[3], CV_RGB(0, 255, 0), 2);
				cvLine(live_show, car[3], car[0], CV_RGB(0, 255, 0), 2);
				cvLine(live_show, car[4], car[5], CV_RGB(255, 0, 0), 2);

				for (int i = 0; i < all_robot_path[i_car_num].size(); i++)
				{
					if (i == 0)
						continue;

					cvLine(live_show, all_robot_path[i_car_num][i - 1], all_robot_path[i_car_num][i], CV_RGB(50, 50, 255), 2);
				}

				for (int i = 0; i < all_robot_path[i_car_num].size(); i++)
				{
					if (i == 0)
						continue;

					if (i % 60 == 0)
					{
						car[0] = cvPoint(all_robot_path[i_car_num][i].x + i_carsize * cos(all_robot_dir[i_car_num][i] + 0.7854), all_robot_path[i_car_num][i].y - i_carsize * sin(all_robot_dir[i_car_num][i] + 0.7854));
						car[1] = cvPoint(all_robot_path[i_car_num][i].x + i_carsize * cos(0.7854 - all_robot_dir[i_car_num][i]), all_robot_path[i_car_num][i].y + i_carsize * sin(0.7854 - all_robot_dir[i_car_num][i]));
						car[2] = cvPoint(all_robot_path[i_car_num][i].x - i_carsize * cos(all_robot_dir[i_car_num][i] + 0.7854), all_robot_path[i_car_num][i].y + i_carsize * sin(all_robot_dir[i_car_num][i] + 0.7854));
						car[3] = cvPoint(all_robot_path[i_car_num][i].x - i_carsize * cos(0.7854 - all_robot_dir[i_car_num][i]), all_robot_path[i_car_num][i].y - i_carsize * sin(0.7854 - all_robot_dir[i_car_num][i]));
						car[4] = cvPoint(all_robot_path[i_car_num][i].x, all_robot_path[i_car_num][i].y);
						car[5] = cvPoint(all_robot_path[i_car_num][i].x + 20 * cos(all_robot_dir[i_car_num][i]), all_robot_path[i_car_num][i].y - 20 * sin(all_robot_dir[i_car_num][i]));


						cvLine(live_show, car[0], car[1], CV_RGB(0, 255, 0), 2);
						cvLine(live_show, car[1], car[2], CV_RGB(0, 255, 0), 2);
						cvLine(live_show, car[2], car[3], CV_RGB(0, 255, 0), 2);
						cvLine(live_show, car[3], car[0], CV_RGB(0, 255, 0), 2);
						cvLine(live_show, car[4], car[5], CV_RGB(255, 0, 0), 2);
					}

				}

				cvPutText(live_show, "Master", cvPoint(i_robot_start_point[i_car_num].x + 10, i_robot_start_point[i_car_num].y - 10), &Font1, Scalar(100, 100, 100));*/
	}
	else
	{
		cvLine(live_show, car[0], car[1], CV_RGB(0, 255, 0), 2);
		cvLine(live_show, car[1], car[2], CV_RGB(0, 255, 0), 2);
		cvLine(live_show, car[2], car[3], CV_RGB(0, 255, 0), 2);
		cvLine(live_show, car[3], car[0], CV_RGB(0, 255, 0), 2);
		cvLine(live_show, car[4], car[5], CV_RGB(255, 150, 150), 2);
		if (detect_area[i_car_num][0].x != 0)
		{
			cvLine(live_show, car[4], detect_area[i_car_num][0], CV_RGB(255, 150, 150), 1);
			cvLine(live_show, car[4], detect_area[i_car_num][1], CV_RGB(255, 150, 150), 1);
			cvLine(live_show, detect_area[i_car_num][0], detect_area[i_car_num][1], CV_RGB(255, 150, 150), 1);
		}

		for (int i = all_robot_path[i_car_num].size(); i > 1; i--)
		{
			if (i == all_robot_path[i_car_num].size() - 500)
				break;
			if (i_car_num == 1)
				cvLine(live_show, all_robot_path[i_car_num][i - 1], all_robot_path[i_car_num][i - 2], CV_RGB(255, 0, 0), 2);
			else if (i_car_num == 2 && total_number > 2)
				cvLine(live_show, all_robot_path[i_car_num][i - 1], all_robot_path[i_car_num][i - 2], CV_RGB(100, 100, 100), 2);
			else if (i_car_num == 3 && total_number > 3)
				cvLine(live_show, all_robot_path[i_car_num][i - 1], all_robot_path[i_car_num][i - 2], CV_RGB(220, 220, 0), 2);
			else if (i_car_num == 4 && total_number > 4)
				cvLine(live_show, all_robot_path[i_car_num][i - 1], all_robot_path[i_car_num][i - 2], CV_RGB(170, 120, 50), 2);
		}
		sprintf_s(num1, "Slave %d", i_car_num);  //機器人編號
		cvPutText(live_show, num1, cvPoint(i_robot_start_point[i_car_num].x + 10, i_robot_start_point[i_car_num].y - 10), &Font1, Scalar(100, 100, 100));
		sprintf_s(num1, "%d", i_arrive_goal_times[i_car_num]); //個別到達次數
		cvPutText(live_show, num1, cvPoint(i_robot_start_point[i_car_num].x + 20, i_robot_start_point[i_car_num].y), &Font2, Scalar(230, 20, 20));

		sprintf_s(num1, "%.4f", show_path_efficiency[i_car_num]);
		cvPutText(live_show, num1, cvPoint(i_robot_start_point[i_car_num].x + 20, i_robot_start_point[i_car_num].y + 30), &Font2, Scalar(230, 20, 20));

	}
	for (int i = 1; i < total_number; i++)
	{
		total_arrive_times += i_arrive_goal_times[i];//計算總共到達次數

		sprintf_s(num1, "%d", car_front_num[i_car_num][i]);
		cvPutText(live_show, num1, cvPoint(i_robot_start_point[i_car_num].x + i * 20, i_robot_start_point[i_car_num].y + 15), &Font2, Scalar(230, 20, 20));


	}

	sprintf_s(num1, "total_arrive_times %d", total_arrive_times);
	cvPutText(live_show, num1, cvPoint(700, 20), &Font2, Scalar(230, 20, 20));
	//	app_pos_m_sim.close();
}

void CCollisionFreepathplanningDlg::car_crash_detector(CvPoint2D64f i_robot_start_point[total_number], int i_carsize, int &io_crash_times, char io_crash_num[total_number][total_number])
{
	for (int i = 1; i < total_number; i++)
	{
		for (int j = i + 1; j < total_number; j++)
		{
			CvPoint2D64f car_temp1 = i_robot_start_point[i], car_temp2 = i_robot_start_point[j];
			double car_dis = sqrt(pow(car_temp1.x - car_temp2.x, 2) + pow(car_temp1.y - car_temp2.y, 2));

			if (car_dis <= i_carsize * 1.9)
			{
				io_crash_num[i][j]++;
				io_crash_num[j][i]++;
				if (io_crash_num[i][j] == 1 || io_crash_num[j][i] == 1)
					io_crash_times++;
			}
			else
			{
				io_crash_num[i][j] = 0;
				io_crash_num[j][i] = 0;
			}
		}
	}
}

void CCollisionFreepathplanningDlg::car_front_detector(CvPoint2D64f i_robot_start_point[total_number], double i_robot_zdir[total_number], double i_detect_angle, int i_detect_radius, CvPoint(&o_detect_area)[total_number][2])
{
	for (int i = 1; i < total_number; i++)
	{
		o_detect_area[i][0] = cvPoint(i_robot_start_point[i].x + i_detect_radius * cos(robot_zdir[i] + i_detect_angle), i_robot_start_point[i].y - i_detect_radius * sin(robot_zdir[i] + i_detect_angle));
		o_detect_area[i][1] = cvPoint(i_robot_start_point[i].x + i_detect_radius * cos(i_detect_angle - robot_zdir[i]), i_robot_start_point[i].y + i_detect_radius * sin(i_detect_angle - robot_zdir[i]));

		CvPoint2D64f car_temp1 = i_robot_start_point[i], car_direction_point = cvPoint2D64f(i_robot_start_point[i].x + 20 * cos(robot_zdir[i]), i_robot_start_point[i].y - 20 * sin(robot_zdir[i]));
		double dx1 = car_direction_point.x - car_temp1.x;
		double dy1 = car_direction_point.y - car_temp1.y;

		for (int j = 1; j < total_number; j++)
		{
			if (i == j)
				continue;

			CvPoint2D64f  car_temp2 = i_robot_start_point[j];
			double car_dis = sqrt(pow(car_temp1.x - car_temp2.x, 2) + pow(car_temp1.y - car_temp2.y, 2));
			double dx2 = car_temp2.x - car_temp1.x;
			double dy2 = car_temp2.y - car_temp1.y;
			double vector_angle = acos((dx1*dx2 + dy1 * dy2) / sqrt((dx1*dx1 + dy1 * dy1)*(dx2*dx2 + dy2 * dy2) + 1e-10));

			if (car_dis < i_detect_radius + 1 && vector_angle < i_detect_angle + 0.01)
				car_front_num[i][j] = true;
			else
				car_front_num[i][j] = false;
		}
	}
}



void CCollisionFreepathplanningDlg::car_path_efficiency(CvPoint2D64f i_robot_start_point, CvPoint2D64f i_robot_end_point, int reentry_tag, int old_reentry_tag, int i_car_num)
{
	double ground_path_length = sqrt(pow(i_robot_start_point.x - i_robot_end_point.x, 2) + pow(i_robot_start_point.y - i_robot_end_point.y, 2)) - 5;
	double real_path_length = 0;
	double temp_path = 0;
	double jump_caculate = 20;
	int i = 0;

	for (i = old_reentry_tag; i < reentry_tag - jump_caculate - 1; i += jump_caculate)
	{
		if (all_robot_path_double[i_car_num][i].x == all_robot_path_double[i_car_num][i + jump_caculate].x && all_robot_path_double[i_car_num][i].y == all_robot_path_double[i_car_num][i + jump_caculate].y)
			continue;

		temp_path = sqrt(pow(all_robot_path_double[i_car_num][i].x - all_robot_path_double[i_car_num][i + jump_caculate].x, 2) + pow(all_robot_path_double[i_car_num][i].y - all_robot_path_double[i_car_num][i + jump_caculate].y, 2)) + temp_path;
	}

	temp_path = sqrt(pow(all_robot_path_double[i_car_num][i].x - all_robot_path_double[i_car_num][reentry_tag - 1].x, 2) + pow(all_robot_path_double[i_car_num][i].y - all_robot_path_double[i_car_num][reentry_tag - 1].y, 2)) + temp_path;

	show_path_efficiency[i_car_num] = ground_path_length / temp_path;
	if (show_path_efficiency[i_car_num] > 1)
		show_path_efficiency[i_car_num] = 1;
}

void CCollisionFreepathplanningDlg::find_path(int x)   // 印出由起點到x點的最短路徑
{
	//	fstream app_truepath_output("路徑輸出.txt", ios::app);

	if (x != parent[x]) // 先把之前的路徑都印出來
		find_path(parent[x]);

	//	app_truepath_output << "x = " << x << " parent[x] = " << parent[x] << endl;
	//	app_truepath_output.close();
	show_path.push_back(x);
}

void CCollisionFreepathplanningDlg::dijkstra(int source, int node_num)
{
	for (int i = 0; i < node_num; i++) visit[i] = false;   // initialize
	for (int i = 0; i < node_num; i++) d[i] = 100000000;

	d[source] = 0;   //令d[a]是起點到a點的最短路徑長度，起點設為零，其他點都是空的
	parent[source] = source;  // 紀錄各個點在最短路徑樹上的父親是誰

	for (int k = 0; k < node_num; k++)
	{
		int a = -1, b = -1, min = 100000000;
		for (int i = 0; i < node_num; i++)
			if (!visit[i] && d[i] < min)
			{
				a = i;  // 記錄這一條邊
				min = d[i];
			}

		if (a == -1) break;     // 起點有連通的最短路徑都已找完
		if (min == 1e9) break;  // 不連通即是最短路徑長度無限長
		visit[a] = true;

		// 以邊ab進行relaxation
		for (b = 0; b < node_num; b++)
			if (!visit[b] && d[a] + w[a][b] < d[b])
			{
				d[b] = d[a] + w[a][b];
				parent[b] = a;
			}
	}
}

void CCollisionFreepathplanningDlg::OnClose()
{
	if (CanExit())
		CDialogEx::OnClose();
}

void CCollisionFreepathplanningDlg::OnOK()
{
	if (CanExit())
		CDialogEx::OnOK();
}

void CCollisionFreepathplanningDlg::OnCancel()
{
	if (CanExit())
		CDialogEx::OnCancel();
}

BOOL CCollisionFreepathplanningDlg::CanExit()
{
	// 如果 Proxy 物件仍在附近，則 Automation 控制器
	// 仍掌控此應用程式。請將對話方塊保持在附近，
	// 但是隱藏其 UI。
	if (m_pAutoProxy != nullptr)
	{
		ShowWindow(SW_HIDE);
		return FALSE;
	}

	return TRUE;
}



void CCollisionFreepathplanningDlg::OnBnClickedOk()
{
	// TODO: 在此加入控制項告知處理常式程式碼
	stop_program = 1;
	CDialogEx::OnOK();
}
