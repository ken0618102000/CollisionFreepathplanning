
// CollisionFree_pathplanningDlg.cpp: 實作檔案
//

#include "stdafx.h"
#include "CollisionFree_pathplanning.h"
#include "CollisionFree_pathplanningDlg.h"
#include "DlgProxy.h"
#include "afxdialogex.h"
#include "fstream"
#include <windows.h>
#include "math.h"
#include "queue"
//#include "CvvImage.cpp"
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
	, check_path_change(FALSE)
	, new_car_wait(FALSE)
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
	DDX_Check(pDX, IDC_CHECK_path_change, check_path_change);
	DDX_Control(pDX, IDC_CHECK_path_change, check_path_change2);
	DDX_Check(pDX, IDC_CHECK_new_car_wait, new_car_wait);
	DDX_Control(pDX, IDC_CHECK_new_car_wait, new_car_wait2);
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
	new_car_wait2.SetCheck(TRUE);

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

int choose_path_num = 0; //同一台機器人中可以選擇要第幾順位之最短路徑
int image_size = 88;
int detecte_dis = 0;
int car_crash_times = 0;
int run_times_counter = 0;
double detect_angle = 0;
double robot_zdir[15] = { 0 };
double show_path_efficiency[total_number] = { 0 };
double show_path_weight[10] = { 0 };
vector <int> show_path[20];
vector <int> path_optimization[total_number];
vector <Point> last_path_optimization[total_number];
vector <Point> path_optimization_real_xy[total_number];
vector <Point> path_optimization_real_xy_for_wait_mechanism[total_number];
vector <Point> all_robot_path[total_number];
vector <Point2d> all_robot_path_double[total_number];
vector <double> all_robot_dir[total_number];
Point2d fix_robot_start_point[total_number] = { 0 }, fix_robot_end_point[total_number] = { 0 };
Point2d delay_own_obstacle[total_number][2] = { 0 };
bool stop_program = false;
bool car_front_num[total_number][total_number] = { false };
char crash_num[total_number][total_number] = { 0 };
voronoi_diagram<double> vd;
rect_type brect_;
vector<point_type> point_data_;
vector<segment_type> segment_data_;

struct Robot_posture_struct
{
	Point2d robot_start_point[total_number] = { 0 },
		robot_end_point[total_number] = { 0 };
	double
		robot_start_direction[total_number] = { 0 },
		robot_end_direction[total_number] = { 0 };

}robot_posture;

struct voro_Point
{
	int a;
	int b;
	voro_Point(int x, int y) : a(x), b(y) {}
};

struct Segment
{
	voro_Point p0;
	voro_Point p1;
	Segment(int x1, int y1, int x2, int y2) : p0(x1, y1), p1(x2, y2) {}
};

namespace boost
{
	namespace polygon
	{

		template <>
		struct geometry_concept<voro_Point>
		{
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
		struct geometry_concept<Segment>
		{
			typedef segment_concept type;
		};

		template <>
		struct segment_traits<Segment>
		{
			typedef int coordinate_type;
			typedef voro_Point point_type;

			static inline point_type get(const Segment& segment, direction_1d dir)
			{
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
	remove("路徑輸出.txt");
	remove("子機器人輸出.txt");
	remove("master_recoder.avi");
	remove("slave_recoder_1.avi");
	remove("slave_recoder_2.avi");
	remove("slave_recoder_3.avi");
	remove("slave_recoder_4.avi");
	remove("slave_recoder_5.avi");
	remove("slave_recoder_6.avi");
	remove("slave_recoder_7.avi");
	remove("slave_recoder_8.avi");
	remove("path_num.txt");

// 	IplConvKernel* pKernel_small = NULL;
// 	pKernel_small = cvCreateStructuringElementEx(11, 11, 5, 5, CV_SHAPE_RECT, NULL);
	Mat pKernel_small = getStructuringElement(MORPH_RECT, Size(11, 11));
	vector<vector<bool>>  sca_image;  //縮圖後二值化的結果
	vector <int> path_optimization;
	vector <Point> save_coner;
	vector <Point> all_point_map;
	vector <Point>  jump_path_optimization, jump_path_optimization_simulation;
	vector <Point2d> all_point_map_original;
	vector <draw_car> car_simulation;
	int corner_count = 0;
	int line_count = 0;  //有多少VD線段
	int new_input_index = 0;
	int path_optimization_size_change = 0;
	bool image_change = 0;
	double Data[3000];
	Point2f* point2 = 0;
	Point2d savepoint1[1500] = { 0.0 }, savepoint2[1500] = { 0.0 };
	Point2d new_savepoint1[1500] = { 0.0 }, new_savepoint2[1500] = { 0.0 };
	LARGE_INTEGER tStart, tEnd, ts;

	Mat read_data = imread("photo/格點化地圖.bmp", 1);
	Mat pGrayImg = Mat::zeros(read_data.size(), CV_8U);
	Mat draw_data = Mat::zeros(read_data.size(), CV_8UC3);
	Mat check_change = Mat::zeros(read_data.size(), CV_8U);
	Mat show_data = Mat::zeros(image_size, image_size, CV_8U); //縮小十倍的矩陣
	Mat read_data_old = Mat::zeros(880, 880, CV_8U);

	m_show2.SetWindowPos(&wndTop, 10, 10, draw_data.cols, draw_data.rows, SWP_SHOWWINDOW);
	char path0[100];
	int photo_conunt = 0;

	int circle_priority[total_target] = { 0, 1, 5, 3, 6, 2, 7, 4, 8 }; //設定朝向角順位
	int circle_priority_T[total_target] = { 0 };
	int circle_priority_opposite[9];

	Point2d default_circle_pos[total_target] = { 0 };
	default_circle_pos[1] = Point2d(120, 340);
	default_circle_pos[2] = Point2d(196, 156);
	default_circle_pos[3] = Point2d(380, 80);
	default_circle_pos[4] = Point2d(564, 156);
	default_circle_pos[5] = Point2d(640, 340);
	default_circle_pos[6] = Point2d(564, 524);
	default_circle_pos[7] = Point2d(380, 600);
	default_circle_pos[8] = Point2d(196, 524);

	double default_circle_dir[total_target] = { 0 };
	default_circle_dir[1] = 0;
	default_circle_dir[2] = -0.25 * CV_PI;
	default_circle_dir[3] = -0.5 * CV_PI;
	default_circle_dir[4] = -0.75 * CV_PI;
	default_circle_dir[5] = -CV_PI;
	default_circle_dir[6] = 0.75 * CV_PI;
	default_circle_dir[7] = 0.5 * CV_PI;
	default_circle_dir[8] = 0.25 * CV_PI;

	for (int i = 1; i < total_target; i++)
	{
		circle_priority_T[circle_priority[i]] = i;
		robot_zdir[circle_priority[i]] = default_circle_dir[i]; //設定朝向角

		if (i < (total_target + 1) / 2) //依照初始點位置自動設定目標點位置
			circle_priority_opposite[circle_priority[i]] = i + (total_target - 1) / 2;
		else
			circle_priority_opposite[circle_priority[i]] = i - (total_target - 1) / 2;
	}

	robot_posture.robot_start_point[1] = default_circle_pos[circle_priority_T[1]];  //An algorithm for multi-robot collision-free (Fig.21測試)
	robot_posture.robot_end_point[1] = default_circle_pos[circle_priority_opposite[1]];

	// 	robot_posture.robot_start_point[1] = Point2d(100, 150);  //An algorithm for multi-robot collision-free (自由測試)
	// 	robot_posture.robot_end_point[1] = Point2d(700, 800);

	if (total_number > 2)
	{
		robot_posture.robot_start_point[2] = default_circle_pos[circle_priority_T[2]];  //An algorithm for multi-robot collision-free (Fig.21測試)
		robot_posture.robot_end_point[2] = default_circle_pos[circle_priority_opposite[2]];

		// 		robot_posture.robot_start_point[2] = Point2d(700, 150);
		// 		robot_posture.robot_end_point[2] = Point2d(100, 800);
	}

	if (total_number > 3)
	{
		robot_posture.robot_start_point[3] = default_circle_pos[circle_priority_T[3]];  //An algorithm for multi-robot collision-free  (Fig.21測試)
		robot_posture.robot_end_point[3] = default_circle_pos[circle_priority_opposite[3]];

		// 		robot_posture.robot_start_point[3] = Point2d(800, 450);
		// 		robot_posture.robot_end_point[3] = Point2d(300, 800);
	}

	if (total_number > 4)
	{
		robot_posture.robot_start_point[4] = default_circle_pos[circle_priority_T[4]];  //An algorithm for multi-robot collision-free (Fig.21測試) 
		robot_posture.robot_end_point[4] = default_circle_pos[circle_priority_opposite[4]];

		// 		robot_posture.robot_start_point[4] = Point2d(700, 700);
		// 		robot_posture.robot_end_point[4] = Point2d(100, 100);

	}

	if (total_number > 5)
	{
		robot_posture.robot_start_point[5] = default_circle_pos[circle_priority_T[5]];   //An algorithm for multi-robot collision-free (Fig.21測試)
		robot_posture.robot_end_point[5] = default_circle_pos[circle_priority_opposite[5]];
	}

	if (total_number > 6)
	{
		robot_posture.robot_start_point[6] = default_circle_pos[circle_priority_T[6]];  //An algorithm for multi-robot collision-free (Fig.21測試)
		robot_posture.robot_end_point[6] = default_circle_pos[circle_priority_opposite[6]];
	}

	if (total_number > 7)
	{
		robot_posture.robot_start_point[7] = default_circle_pos[circle_priority_T[7]];   //An algorithm for multi-robot collision-free (Fig.21測試)
		robot_posture.robot_end_point[7] = default_circle_pos[circle_priority_opposite[7]];
	}
	if (total_number > 8)
	{
		robot_posture.robot_start_point[8] = default_circle_pos[circle_priority_T[8]];   //An algorithm for multi-robot collision-free (Fig.21測試)
		robot_posture.robot_end_point[8] = default_circle_pos[circle_priority_opposite[8]];
	}
	if (total_number > 9)
	{
		robot_posture.robot_start_point[9].x = 160;  //路徑起始與終點，請參照圖片給定
		robot_posture.robot_start_point[9].y = 230;
		robot_posture.robot_end_point[9].x = 600;
		robot_posture.robot_end_point[9].y = 480;
	}
	if (total_number > 10)
	{
		robot_posture.robot_start_point[10].x = 600;  //路徑起始與終點，請參照圖片給定
		robot_posture.robot_start_point[10].y = 480;
		robot_posture.robot_end_point[10].x = 160;
		robot_posture.robot_end_point[10].y = 230;
	}
	if (total_number > 11)
	{
		robot_posture.robot_start_point[11].x = 250;  //路徑起始與終點，請參照圖片給定
		robot_posture.robot_start_point[11].y = 120;
		robot_posture.robot_end_point[11].x = 510;
		robot_posture.robot_end_point[11].y = 560;
	}
	if (total_number > 12)
	{
		robot_posture.robot_start_point[12].x = 510;  //路徑起始與終點，請參照圖片給定
		robot_posture.robot_start_point[12].y = 560;
		robot_posture.robot_end_point[12].x = 250;
		robot_posture.robot_end_point[12].y = 120;
	}
	if (total_number > 13)
	{
		robot_posture.robot_start_point[13].x = 800;  //路徑起始與終點，請參照圖片給定
		robot_posture.robot_start_point[13].y = 250;
		robot_posture.robot_end_point[13].x = 100;
		robot_posture.robot_end_point[13].y = 650;
	}
	memcpy(fix_robot_start_point, robot_posture.robot_start_point, sizeof(fix_robot_start_point));
	memcpy(fix_robot_end_point, robot_posture.robot_end_point, sizeof(fix_robot_end_point));

	photo_conunt = 1;
	while (true)
	{
		QueryPerformanceFrequency(&ts);
		QueryPerformanceCounter(&tStart);

		sprintf_s(path0, "photo/空白.png", photo_conunt);
		fstream in_image0(path0, ios::in);

		if (!in_image0)
			break;

		read_data = imread(path0, 0);

		if (1)
		{
			//-------------------------清除數據------------------------------
			save_coner.clear();
			sca_image.clear();
			path_optimization.clear();
			all_point_map.clear();
			all_point_map_original.clear();
			for (int i = 0; i < 20; i++)
				show_path[i].clear();
			memset(show_path_weight, 0, sizeof(show_path_weight));
			draw_data.zeros(draw_data.rows, draw_data.cols, CV_8UC3);
//			memset((unsigned char*)draw_data->imageData, 0, draw_data->imageSize);
			//---------------------------------------------------------------------
			resize(read_data, pGrayImg, Size(read_data.rows, read_data.cols), 0, 0, INTER_NEAREST);//讀黑白影像用的
	//		cvCvtColor(read_data, pGrayImg, CV_RGB2GRAY);  //讀彩色影像用的
			erode(pGrayImg, pGrayImg, pKernel_small, Point(-1, -1), 3);  //侵蝕的相反(因為是白底)
// 			cvDilate(pGrayImg, pGrayImg, pKernel_small, 2);  //膨脹的相反
			cvtColor(pGrayImg, draw_data, COLOR_GRAY2RGB);
			//			cvSaveImage("給連通物件用的.bmp", pGrayImg);

						//數值要依據縮小倍率與格點pixel數決定


			resize(pGrayImg, show_data, Size(show_data.rows, show_data.cols), 0, 0, INTER_NEAREST);


			//輸入圖片，輸出二值資料
			binarization(show_data, sca_image);
			//輸入二值資料，輸出角點
			find_coner(sca_image, save_coner, 4);
			//將角點轉換為準備要丟入Voronoi運算的格式
			trans2Voronoi(sca_image, save_coner, Data, 8);
			//計算狹義Voronoi，輸入角點資料與邊界，輸出兩個矩陣 

			//-------------------------------------------繪圖---------------------------------------

			MultiRobot_Path_simulation(show_data, jump_path_optimization, sca_image, save_coner, robot_posture.robot_start_point, robot_posture.robot_end_point, 2, jump_path_optimization_simulation, car_simulation, draw_data);  //開始模擬
			jump_path_optimization.clear();

			QueryPerformanceCounter(&tEnd);

			m_total_time = 1000 / ((tEnd.QuadPart - tStart.QuadPart) * 1000 / (double)(ts.QuadPart));
			m_coner_count = Data[0] + 1; //顯示角點數量
			UpdateData(FALSE);
		}


		show_data.zeros(show_data.rows, show_data.cols, CV_8UC3);
// 		memset((unsigned char*)show_data->imageData, 0, show_data->imageSize);
// 		cvReleaseImage(&read_data);


		break;
		//		photo_conunt++;
	}


	//-------------------------清除數據------------------------------
	save_coner.clear();
	sca_image.clear();
	path_optimization.clear();
	all_point_map.clear();
	all_point_map_original.clear();
	for (int i = 0; i < 20; i++)
		show_path[i].clear();
	memset(show_path_weight, 0, sizeof(show_path_weight));
//	memset((unsigned char*)draw_data->imageData, 0, draw_data->imageSize);
	point_data_.clear();
	segment_data_.clear();
	vd.clear();

	//---------------------------------------------------------------------

	CDialogEx::OnOK();
}


void CCollisionFreepathplanningDlg::MultiRobot_Path_simulation(Mat i_draw_data, vector <Point> i_host_path, vector<vector<bool>>  i_sca_image, vector <Point> i_save_coner, Point2d i_robot_start_point[total_number], Point2d i_robot_end_point[total_number], int i_car_density, vector <Point>& o_sim_path, vector <draw_car>& o_sim_car, Mat& offline_show)
{
	VideoWriter master_recoder;
	char master_recoder_name[] = "master_recoder.avi";
	int arrive_times_set = 25;
	int FPS = 25;
	Size AviSize = Size(image_size * 10, image_size * 10);
	int AviColor = 1;
	int Reentry_tag[total_number] = { 0 }, old_Reentry_tag[total_number] = { 0 };
	//	master_recoder = cvCreateVideoWriter(master_recoder_name, CV_FOURCC('D', 'I', 'V', 'X'), FPS, AviSize, AviColor);  //20180926拿掉Master



	VideoWriter slave_recoder[9];
	char slave_recoder_name_1[] = "slave_recoder_1.avi";
	slave_recoder[1] = VideoWriter(slave_recoder_name_1, VideoWriter::fourcc('D', 'I', 'V', 'X'), FPS, AviSize, true);

	char slave_recoder_name_2[] = "slave_recoder_2.avi";
	slave_recoder[2] = VideoWriter(slave_recoder_name_2, VideoWriter::fourcc('D', 'I', 'V', 'X'), FPS, AviSize, true);

	char slave_recoder_name_3[] = "slave_recoder_3.avi";
	slave_recoder[3] = VideoWriter(slave_recoder_name_3, VideoWriter::fourcc('D', 'I', 'V', 'X'), FPS, AviSize, true);

	char slave_recoder_name_4[] = "slave_recoder_4.avi";
	slave_recoder[4] = VideoWriter(slave_recoder_name_4, VideoWriter::fourcc('D', 'I', 'V', 'X'), FPS, AviSize, true);

	char slave_recoder_name_5[] = "slave_recoder_5.avi";
	slave_recoder[5] = VideoWriter(slave_recoder_name_4, VideoWriter::fourcc('D', 'I', 'V', 'X'), FPS, AviSize, true);

	char slave_recoder_name_6[] = "slave_recoder_6.avi";
	slave_recoder[6] = VideoWriter(slave_recoder_name_4, VideoWriter::fourcc('D', 'I', 'V', 'X'), FPS, AviSize, true);

	char slave_recoder_name_7[] = "slave_recoder_7.avi";
	slave_recoder[7] = VideoWriter(slave_recoder_name_4, VideoWriter::fourcc('D', 'I', 'V', 'X'), FPS, AviSize, true);

	char slave_recoder_name_8[] = "slave_recoder_8.avi";
	slave_recoder[8] = VideoWriter(slave_recoder_name_4, VideoWriter::fourcc('D', 'I', 'V', 'X'), FPS, AviSize, true);
	//VideoWriter::fourcc

	double sampleTime = 38;
	double pixel2cm = 100;
	double rho_dot, alpha_dot, beta_dot;
	double x_here = 0, y_here = 0, zdir_here = 0, theta_here;
	double u1, u2;
	double target_pos_sim[3] = { 0 }, car_x_sim, car_y_sim, car_zdir_sim, phi_sim, rho_sim, theta_sim, alpha_sim, beta_sim;
	double vr = 0, vl = 0;
	Point erase_path = { 0 };
	vector <double> x_save, y_save;
	vector <draw_car>local_draw_car;
	int jump_draw = 0, state;
	//Point draw_car[6];
	int carsize = 20;
	int initial_scale = 100;
	int arrive_goal_times[total_number] = { 0 };
	int arrive_time_cost[total_number] = { 0 };

	// 	cvSetZero(draw_data);
	Point draw_oringin[2];

	remove("control_pos_m_sim.txt");
	remove("car_path_efficiency.txt");
	fstream app_car_path_efficiency("car_path_efficiency.txt", ios::app);



	while (true)
	{
		run_times_counter++; // 計數跑的次數
		if (stop_program == 1)
			break;

		int all_robot_arrive_ones = 0;
		for (int robot_number = 1; robot_number < total_number; robot_number++)
		{
			if (arrive_goal_times[robot_number] == 1) //所有機器人走到目標點即停止
			{
				all_robot_arrive_ones++;
				continue;
			}

			if (abs(i_robot_start_point[robot_number].x - i_robot_end_point[robot_number].x) > 2 || abs(i_robot_start_point[robot_number].y - i_robot_end_point[robot_number].y) > 2)
			{
				Servant_Path(robot_number, i_draw_data, i_sca_image, i_save_coner, i_robot_start_point, i_robot_end_point, o_sim_path, carsize, slave_recoder[robot_number], arrive_goal_times);
				arrive_time_cost[robot_number]++;
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
				Servant_Path(robot_number, i_draw_data, i_sca_image, i_save_coner, i_robot_start_point, i_robot_end_point, o_sim_path, carsize, slave_recoder[robot_number], arrive_goal_times);

				app_car_path_efficiency << robot_number << " " << arrive_goal_times[robot_number] << " 路徑效率= " << show_path_efficiency[robot_number] << " 花費時間= " << arrive_time_cost[robot_number] << endl;
			}
		}

		for (int i = 0; i < total_number; i++)
		{
			path_optimization[i].clear();
			path_optimization_real_xy[i].clear();
			path_optimization_real_xy_for_wait_mechanism[i].clear();
		}

		if (arrive_goal_times[1] > arrive_times_set&& arrive_goal_times[2] > arrive_times_set&& arrive_goal_times[3] > arrive_times_set&& arrive_goal_times[4] > arrive_times_set || all_robot_arrive_ones == total_number - 1) //all_robot_arrive_ones == total_number - 1 (跑完一次自動關閉程式)
			break;

	}

	app_car_path_efficiency.close();

	o_sim_car = local_draw_car;
	local_draw_car.clear();
	return;

}

void CCollisionFreepathplanningDlg::Servant_Path(int i_robot_num, Mat i_pGrayImg, vector<vector<bool>> i_sca_image, vector <Point> i_save_coner, Point2d i_robot_start_point[total_number], Point2d i_robot_end_point[total_number], vector <Point> i_sim_path, int i_carsize, VideoWriter i_slave_recoder, int i_arrive_goal_times[total_number])
{
	int delay_update_parameter = 4;  //更新自身voronoi coner的頻率。數值為機器人長度的倍數，1及代表1倍
	int detect_dis_parameter = 2;   //偵測前方機器人的距離參數調整。數值為機器人長度的倍數，1及代表延伸機器人大小之1倍
	double detect_angle_parameter = 90;  //偵測前方機器人的角度範圍。輸入數值為角度

	int corner_count = 0;
	int line_count = 0;  //有多少VD線段
	int new_input_index = 0;
	int path_optimization_size_change = 0;
	int robot_Dilate_size = (i_carsize / 10) + 2;  //看高優先權的機器人要膨脹到多大
	char num1[50];
	double Data[3000] = { 0 };
	Point2f* point2 = 0;
	Point2d savepoint1[1500] = { 0.0 }, savepoint2[1500] = { 0.0 };
	Point2d new_savepoint1[1500] = { 0.0 }, new_savepoint2[1500] = { 0.0 };
	vector <Point> all_point_map;
	vector <Point>  jump_path_optimization, jump_path_optimization_simulation;
	vector <Point> save_coner;
	vector <Point2d> all_point_map_original;

	Mat draw_data6 = Mat(image_size, image_size, CV_8UC3);
	Mat draw_data5 = Mat(image_size * 5, image_size * 5, CV_8UC3);
	Mat draw_data4 = Mat(image_size, image_size, CV_8UC3);
	Mat draw_data3 = i_pGrayImg.clone();
	Mat draw_data2 = Mat(image_size * 10, image_size * 10, CV_8UC3);
	Mat draw_data1 = Mat(image_size, image_size, CV_8U);

	//	draw_data3 = cvCreateImage(cvGetSize(i_pGrayImg), 8, 3);

	Point host_position[50];
	Point other_position;

	//----------------------------清除數據----------------------------------------------------------------------
	for (int i = 0; i < 20; i++)
		show_path[i].clear();
	point_data_.clear();
	vd.clear();
	//-------------------------------------------------------------------------------------------------------------
	delay_own_obstacle[i_robot_num][1] = i_robot_start_point[i_robot_num];
	double own_dis = sqrt(pow(delay_own_obstacle[i_robot_num][1].x - delay_own_obstacle[i_robot_num][0].x, 2) + pow(delay_own_obstacle[i_robot_num][1].y - delay_own_obstacle[i_robot_num][0].y, 2));

	detecte_dis = detect_dis_parameter;
	detect_angle = detect_angle_parameter * CV_PI / 360; //偵測前方機器人的角度範圍。輸入數值為角度
	//新增延遲更新系統，當自己的coner離自己一段距離後才更新到當下位置
	if (own_dis > delay_update_parameter* i_carsize)  //數值為機器人長度的倍數，1及代表1倍
		delay_own_obstacle[i_robot_num][0] = i_robot_start_point[i_robot_num];

	for (int i = i_robot_num; i > 0; i--)  //新增主從機器人之障礙格點
	{
		if (i_robot_start_point[i].x == 0 && i_robot_start_point[i].y == 0 || i == i_robot_num)
			continue;

		other_position.x = i_robot_start_point[i].x / 10;  //看是要用即時修正還是要延遲修正則用不同的公式
		other_position.y = i_robot_start_point[i].y / 10;  //delay_own_obstacle[i][0]延遲更新
																					// i_robot_start_point[i]即時更新
		save_coner.push_back(other_position);

		for (int j = -robot_Dilate_size; j <= robot_Dilate_size; j++)  //看高優先權的機器人要膨脹到多大
		{
			for (int k = -robot_Dilate_size; k <= robot_Dilate_size; k++)
			{

				i_sca_image[other_position.y + j][other_position.x + k] = 1;
				draw_data3.at<uchar>(other_position.y + j, other_position.x + k) = 0;
			}
		}

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
	draw_data3.at<uchar>(other_position.y, other_position.x) = 0;
	//------------------------------------------------------

	//將角點轉換為準備要丟入Voronoi運算的格式
	trans2Voronoi(i_sca_image, save_coner, Data, 8);
	//計算狹義Voronoi，輸入角點資料與邊界，輸出兩個矩陣
//	Voronoi_calculate(Data, image_size, image_size, savepoint1, savepoint2, line_count);
	//計算狹義boost_Voronoi，輸入角點資料與邊界，輸出兩個矩陣，新版使用boost c++
	boost_Voronoi_calculate(image_size, image_size, savepoint1, savepoint2, line_count);
	//計算廣義Voronoi，待改
	Generalized_Voronoi(i_sca_image, savepoint1, savepoint2, line_count, new_input_index, new_savepoint1, new_savepoint2);
	//VD點會破碎，將其重新聚合
//	Match_point(line_count, new_input_index, new_savepoint1, new_savepoint2, 1);
	//Dijkstra路徑搜尋，輸入點連接資訊跟數量
	Dijkstra_path_planning(i_robot_num, i_robot_start_point, i_robot_end_point, new_savepoint1, new_savepoint2, new_input_index, all_point_map, all_point_map_original);

	//------------------------------------------------------
	save_coner.pop_back();  //再刪除自身障礙物來連線
	i_sca_image[other_position.y][other_position.x] = 0;
	draw_data3.at<uchar>(other_position.y, other_position.x) = 1;
	//------------------------------------------------------


	//路徑優化，輸入二值資訊與原本路徑
	Path_Optimization_20181127(i_sca_image, all_point_map_original, i_robot_num, path_optimization[i_robot_num]); //新版捷徑搜尋，可以跳過障礙物之後繼續搜尋路徑


	Point final_goal_point;
	final_goal_point.x = fix_robot_end_point[i_robot_num].x;
	final_goal_point.y = fix_robot_end_point[i_robot_num].y;
	all_point_map.push_back(final_goal_point);

	if (check_path_change2.GetCheck())  //(paper 2)手動切換選擇路徑
	{
		//低優先權的機器人要檢查在未來有沒有機會阻擋到高優先權的機器人，如果有則規劃另一條路線。
		Check_Thepath_Collision(i_sca_image, all_point_map_original, i_robot_num, all_point_map, path_optimization);
	}

  //以座標的方式儲存最新path_optimization後的結果
	last_path_optimization[i_robot_num].clear();
	for (unsigned int path_opt = 0; path_opt < path_optimization[i_robot_num].size(); path_opt++) //路徑優化顯示
	{
		last_path_optimization[i_robot_num].push_back(Point(all_point_map[path_optimization[i_robot_num][path_opt]].x, all_point_map[path_optimization[i_robot_num][path_opt]].y));
	}

	Mat pKernel_small = getStructuringElement(MORPH_RECT, Size(3, 3));

	dilate(draw_data3, draw_data1, pKernel_small, Point(-1, -1), 2);
	//	cvSaveImage("first/output.png", draw_data1);
	cvtColor(draw_data1, draw_data4, COLOR_GRAY2RGB);
	cvtColor(draw_data3, draw_data6, COLOR_GRAY2RGB);
	addWeighted(draw_data4, 0.5, draw_data6, 0.5, 0, draw_data4);  //疊加膨脹後的障礙物，可視化用
	resize(draw_data4, draw_data2, draw_data2.size(), 0, 0, INTER_NEAREST);
	//	cvSaveImage("first/output1.png", draw_data2);
			//-------------------------------------------繪圖---------------------------------------


// 	for (int i = 0; i < save_coner.size(); i++)  //障礙物角點圖
// 	{
// 		line(draw_data2, Point(save_coner[i].x * 10, save_coner[i].y * 10), Point(save_coner[i].x * 10, save_coner[i].y * 10), CV_RGB(255, 0, 0), 12);
// 	}
	for (int i = 0; i < point_data_.size(); i++) //全部角點圖
	{
		line(draw_data2, Point(point_data_[i].x() * 10, point_data_[i].y() * 10), Point(point_data_[i].x() * 10, point_data_[i].y() * 10), CV_RGB(255, 0, 0), 12);
	}
	//imwrite("first/output2.png", draw_data2);
	// 	imwrite("first/output2.png", draw_data2);
	for (int i = 0; i < line_count; i++)   //VD圖
	{
		line(draw_data2, Point(savepoint1[i].x * 10, savepoint1[i].y * 10), Point(savepoint2[i].x * 10, savepoint2[i].y * 10), CV_RGB(0, 0, 255), 1);
	}
	for (int i = 0; i < new_input_index; i++)  //GVD圖
	{
		line(draw_data2, Point(new_savepoint1[i].x * 10, new_savepoint1[i].y * 10), Point(new_savepoint2[i].x * 10, new_savepoint2[i].y * 10), CV_RGB(0, 0, 255), 2);
	}
	// 	imwrite("first/output3.png", draw_data2);
	for (int i = 0; show_path[i].size() != 0; i++) //畫出所有路徑圖
	{
		for (unsigned int path_index = 0; path_index < show_path[i].size() - 1; path_index++)
		{
			line(draw_data2, Point(all_point_map[show_path[i][path_index]].x, all_point_map[show_path[i][path_index]].y), Point(all_point_map[show_path[i][path_index + 1]].x, all_point_map[show_path[i][path_index + 1]].y), CV_RGB(255 - i*20, 50, 5+i*20), 2);
		}
//		break;  //要畫出全部路徑圖就註解掉此行
	}
	// 	imwrite("first/output4.png", draw_data2); 
	for (int i = 0; i < all_point_map.size(); i++)  //顯示所有VD點編號
	{
		sprintf_s(num1, "%d", i);
		putText(draw_data2, num1, Point(all_point_map[i].x + 5, all_point_map[i].y), FONT_HERSHEY_PLAIN, 1, Scalar(170, 170, 170), 1, LINE_AA);
	}

	for (unsigned int path_opt = 0; path_opt < path_optimization[i_robot_num].size(); path_opt++) //路徑優化顯示
	{
		jump_path_optimization.push_back(all_point_map[path_optimization[i_robot_num][path_opt]]);

		if (path_opt == path_optimization[i_robot_num].size() - 1)
			break;
		else
		{
			line(draw_data2, Point(all_point_map[path_optimization[i_robot_num][path_opt]].x, all_point_map[path_optimization[i_robot_num][path_opt]].y), Point(all_point_map[path_optimization[i_robot_num][path_opt + 1]].x, all_point_map[path_optimization[i_robot_num][path_opt + 1]].y), CV_RGB(0, 200, 0), 3);
		}
	}
	//	imwrite("first/output6.png", draw_data2);
	if (jump_path_optimization.size() < 2)
	{
		jump_path_optimization.push_back(Point(i_robot_end_point[i_robot_num].x, i_robot_end_point[i_robot_num].y));
		// 		cvReleaseImage(&draw_data2);  //原為到達點之後就停止，已改為無限循環
		// 		cvReleaseImage(&draw_data1);
		// 		cvReleaseImage(&draw_data3);
		// 		return;
	}

	Point2d temp_xy_double;
	Point temp_xy;

	if (!car_wait_mechanism(i_robot_start_point, robot_zdir, i_robot_num, all_point_map, last_path_optimization) || !new_car_wait2.GetCheck())  //(paper 2)如果沒有不用等待的話就進入算下一步
		ServantRobot_Path_simulation(jump_path_optimization, i_robot_num, temp_xy_double, i_robot_start_point[i_robot_num], robot_zdir[i_robot_num]);
	else
	{
		temp_xy_double = i_robot_start_point[i_robot_num];
		robot_zdir[i_robot_num] = all_robot_dir[i_robot_num][all_robot_dir[i_robot_num].size() - 1]; //如果要等待的話就直接將上一次位置的數據放入
	}

//	ServantRobot_Path_simulation(jump_path_optimization, i_robot_num, temp_xy_double, i_robot_start_point[i_robot_num], robot_zdir[i_robot_num]);

	i_robot_start_point[i_robot_num] = temp_xy_double;
	temp_xy.x = i_robot_start_point[i_robot_num].x;
	temp_xy.y = i_robot_start_point[i_robot_num].y;
	all_robot_path[i_robot_num].push_back(temp_xy);
	all_robot_dir[i_robot_num].push_back(robot_zdir[i_robot_num]);
	all_robot_path_double[i_robot_num].push_back(temp_xy_double);

	for (int i = 0; i < total_number; i++)
		simulation_car(draw_data2, i_robot_start_point, robot_zdir, i, i_carsize, i_arrive_goal_times);

	//	imwrite("first/output5.png", draw_data2);

	car_crash_detector(i_robot_start_point, i_carsize, car_crash_times, crash_num);  //偵測有無碰撞
	sprintf_s(num1, "Crash times %d", car_crash_times);
	putText(draw_data2, num1, Point(30, 30), FONT_HERSHEY_PLAIN, 1, Scalar(50, 50, 255), 1, LINE_AA);

	if (check_path_change2.GetCheck())  
	{
		sprintf_s(num1, "paper %d", 2); //顯示新版資訊
		putText(draw_data2, num1, Point(700, 20), 1, 1, Scalar(230, 20, 20), 1, LINE_8);
	}

	if (i_robot_num == 1)
	{
		imshow("draw_data", draw_data2); // 動態顯示模擬
		waitKey(5);
		show_path;
		show_path_weight;

//		i_slave_recoder.write(draw_data2);
	}

//	resize(draw_data2, draw_data5, Size(draw_data5.rows, draw_data5.cols), 0, 0, INTER_NEAREST);
	i_slave_recoder.write(draw_data2);

}

void CCollisionFreepathplanningDlg::ServantRobot_Path_simulation(vector <Point> i_Servant_path, int i_robot_num, Point2d& o_ServantRobot_pos, Point2d i_robot_start_point, double& io_zdir)
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
	//Point draw_car[6];
	int initial_scale = 100;

	// 	cvSetZero(draw_data);

	// 	remove("control_pos_m_sim.txt");
	// 	fstream app_pos_m_sim("control_pos_m_sim.txt", ios::app);

	// 	draw_oringin[0] = Point(orgin.x, orgin.y);
	// 	draw_oringin[1] = Point(orgin.x, 700 - orgin.y);

	vector <Point>  jump_path_optimization_copy_sim;
	jump_path_optimization_copy_sim.assign(i_Servant_path.begin(), i_Servant_path.end());


	double x_transpos = jump_path_optimization_copy_sim[1].x;
	double y_transpos = jump_path_optimization_copy_sim[1].y;


	double distant_check = sqrt((i_robot_start_point.x - jump_path_optimization_copy_sim[1].x - target_pos_sim[0]) * (i_robot_start_point.x - jump_path_optimization_copy_sim[1].x - target_pos_sim[0]) + (-i_robot_start_point.y + jump_path_optimization_copy_sim[1].y - target_pos_sim[1]) * (-i_robot_start_point.y + jump_path_optimization_copy_sim[1].y - target_pos_sim[1]));



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
	car_zdir_sim = io_zdir;/*+ 0.148+ (double)rand() / (RAND_MAX + 1.0) * 6*/;   //要更改初始方向角的話在這裡加   /*Camera[5] +*/

	if (distant_check > 50)
	{
		target_pos_sim[0] = 5 * car_x_sim / 6;
		target_pos_sim[1] = 5 * car_y_sim / 6;
	}

	phi_sim = car_zdir_sim;
	if (phi_sim > CV_PI)		phi_sim = -2 * CV_PI + phi_sim;
	if (phi_sim < -CV_PI)		phi_sim = 2 * CV_PI + phi_sim;

	rho_sim = sqrt((car_x_sim - target_pos_sim[0]) * (car_x_sim - target_pos_sim[0]) + (car_y_sim - target_pos_sim[1]) * (car_y_sim - target_pos_sim[1])) * pixel2cm / 100;
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

	if (!new_car_wait2.GetCheck())
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
	alpha_dot = (sin(alpha_sim) / rho_sim) * u1 - u2;
	beta_dot = -(sin(alpha_sim) / rho_sim) * u1;

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

void CCollisionFreepathplanningDlg::Path_Optimization_20181127(vector<vector<bool>> i_sca_image, vector<Point2d> i_all_point_map_original, int i_robot_num, vector<int>& o_path_optimization)
{

	int path_inside, break_index, line_distant, crash_index;
	unsigned int cheak_x, cheak_y, cheak_collision_index = 0, final_noncollision_index = 0;
	Point2d cutout_point[500];
	Point2d temp;
	vector <bool> allpath_crash_index;

	temp.x = fix_robot_end_point[i_robot_num].x / 10;
	temp.y = fix_robot_end_point[i_robot_num].y / 10;

	o_path_optimization.push_back(show_path[choose_path_num][0]);
	show_path[choose_path_num].push_back(i_all_point_map_original.size());
	i_all_point_map_original.push_back(temp);

	for (unsigned int i = 0; i < show_path[choose_path_num].size() - 1;)
	{
		crash_index = 0;
		for (unsigned int j = i + 1; j < show_path[choose_path_num].size(); j++)
		{
			line_distant = hypot(i_all_point_map_original[show_path[choose_path_num][i]].x - i_all_point_map_original[show_path[choose_path_num][j]].x, i_all_point_map_original[show_path[choose_path_num][i]].y - i_all_point_map_original[show_path[choose_path_num][j]].y);
//			line_distant = sqrt(pow(i_all_point_map_original[show_path[choose_path_num][i]].x - i_all_point_map_original[show_path[choose_path_num][j]].x, 2) + pow(i_all_point_map_original[show_path[choose_path_num][i]].y - i_all_point_map_original[show_path[choose_path_num][j]].y, 2));

			int cutout = 0;
			for (cutout = 0; cutout < line_distant; cutout++)
			{
				//計算斜線上的點
				if (i_all_point_map_original[show_path[choose_path_num][i]].x >= i_all_point_map_original[show_path[choose_path_num][j]].x && i_all_point_map_original[show_path[choose_path_num][i]].y >= i_all_point_map_original[show_path[choose_path_num][j]].y)
				{
					cutout_point[cutout].x = i_all_point_map_original[show_path[choose_path_num][j]].x + (sqrt(pow((i_all_point_map_original[show_path[choose_path_num][i]].x - i_all_point_map_original[show_path[choose_path_num][j]].x), 2)) / (double)line_distant) * (double)cutout;
					cutout_point[cutout].y = i_all_point_map_original[show_path[choose_path_num][j]].y + (sqrt(pow((i_all_point_map_original[show_path[choose_path_num][i]].y - i_all_point_map_original[show_path[choose_path_num][j]].y), 2)) / (double)line_distant) * (double)cutout;
				}

				if (i_all_point_map_original[show_path[choose_path_num][i]].x >= i_all_point_map_original[show_path[choose_path_num][j]].x && i_all_point_map_original[show_path[choose_path_num][j]].y >= i_all_point_map_original[show_path[choose_path_num][i]].y)
				{
					cutout_point[cutout].x = i_all_point_map_original[show_path[choose_path_num][j]].x + (sqrt(pow((i_all_point_map_original[show_path[choose_path_num][i]].x - i_all_point_map_original[show_path[choose_path_num][j]].x), 2)) / (double)line_distant) * (double)cutout;
					cutout_point[cutout].y = i_all_point_map_original[show_path[choose_path_num][j]].y - (sqrt(pow((i_all_point_map_original[show_path[choose_path_num][i]].y - i_all_point_map_original[show_path[choose_path_num][j]].y), 2)) / (double)line_distant) * (double)cutout;
				}

				if (i_all_point_map_original[show_path[choose_path_num][j]].x >= i_all_point_map_original[show_path[choose_path_num][i]].x && i_all_point_map_original[show_path[choose_path_num][i]].y >= i_all_point_map_original[show_path[choose_path_num][j]].y)
				{
					cutout_point[cutout].x = i_all_point_map_original[show_path[choose_path_num][j]].x - (sqrt(pow((i_all_point_map_original[show_path[choose_path_num][i]].x - i_all_point_map_original[show_path[choose_path_num][j]].x), 2)) / (double)line_distant) * (double)cutout;
					cutout_point[cutout].y = i_all_point_map_original[show_path[choose_path_num][j]].y + (sqrt(pow((i_all_point_map_original[show_path[choose_path_num][i]].y - i_all_point_map_original[show_path[choose_path_num][j]].y), 2)) / (double)line_distant) * (double)cutout;
				}

				if (i_all_point_map_original[show_path[choose_path_num][j]].x >= i_all_point_map_original[show_path[choose_path_num][i]].x && i_all_point_map_original[show_path[choose_path_num][j]].y >= i_all_point_map_original[show_path[choose_path_num][i]].y)
				{
					cutout_point[cutout].x = i_all_point_map_original[show_path[choose_path_num][i]].x + (sqrt(pow((i_all_point_map_original[show_path[choose_path_num][i]].x - i_all_point_map_original[show_path[choose_path_num][j]].x), 2)) / (double)line_distant) * (double)cutout;
					cutout_point[cutout].y = i_all_point_map_original[show_path[choose_path_num][i]].y + (sqrt(pow((i_all_point_map_original[show_path[choose_path_num][i]].y - i_all_point_map_original[show_path[choose_path_num][j]].y), 2)) / (double)line_distant) * (double)cutout;
				}
			}

			for (int cheak_pixel = 0; cheak_pixel < cutout; cheak_pixel++) //檢查兩點間的路徑上是否有障礙物並建立allpath_crash_index表
			{
				cheak_x = round(cutout_point[cheak_pixel].x);
				cheak_y = round(cutout_point[cheak_pixel].y);

				if ((cheak_y > i_sca_image.size() - 2) || (cheak_x > i_sca_image.size() - 2) || (cheak_y < 1 || (cheak_x < 1)))
					continue;

				if (i_sca_image[cheak_y][cheak_x] == 1 /*||
					i_sca_image[cheak_y + 1][cheak_x] == 1 ||
					i_sca_image[cheak_y - 1][cheak_x] == 1 ||
					i_sca_image[cheak_y][cheak_x + 1] == 1 ||
					i_sca_image[cheak_y][cheak_x - 1] == 1 ||
					i_sca_image[cheak_y - 1][cheak_x - 1] == 1 ||
					i_sca_image[cheak_y + 1][cheak_x + 1] == 1 ||
					i_sca_image[cheak_y - 1][cheak_x + 1] == 1 ||
					i_sca_image[cheak_y + 1][cheak_x - 1] == 1*/
					)
				{
					allpath_crash_index.push_back(1);  //如果有障礙物為1
					crash_index = 1;
					break;
				}
				else
				{
					crash_index = 0;
				}
			}

			if (crash_index == 0)
				allpath_crash_index.push_back(0); //如果有障礙物為0

		}

		for (int k = allpath_crash_index.size() - 1; k >= 0; k--)  //從末端逐一檢查路徑
		{
			if (allpath_crash_index[k] == 0)  //如果與最後一點有路徑可走就直接當作捷徑並結束整個運算
			{
				cheak_collision_index = k + 1;
				break;
			}
			else if (k == 0)  //假如都沒有捷徑的話就檢查下一個路徑點跟剩下的路徑點
			{
				cheak_collision_index = 2;
			}
			else if (k == i)  //逆向檢查、直到檢查到無碰撞的那條捷徑
			{
				cheak_collision_index = i + 1;
				break;
			}
			else
			{
				allpath_crash_index.pop_back();  //從尾部剔除allpath_crash_index裡會碰撞的路徑
			}
		}
		o_path_optimization.push_back(show_path[choose_path_num][cheak_collision_index]); //將捷徑傳回o_path_optimization

		i = cheak_collision_index;  //檢查剩餘路段是否有捷徑
		cheak_collision_index = 0;
	}
}


void CCollisionFreepathplanningDlg::Check_Thepath_Collision(vector<vector<bool>> i_sca_image, vector <Point2d> i_all_point_map_original, int i_robot_num, vector <Point> i_all_point_map, vector <int> i_path_optimization[total_number])
{
	//(paper 2)
	//準備要將show_path裡的點的座標抓出來
	int crash_with_num = 0, contact_distan_with_high_priority_robot[total_number] = { 0 }, contact_distan_with_low_priority_robot[total_number] = { 0 };
	bool is_the_line_intersect[total_number] = { 0 };
	Point current_robot_line[2], high_priority_robot_line[2], have_intersection[total_number];

	for (int i = 0; i < i_path_optimization[i_robot_num].size(); i++)
		path_optimization_real_xy[i_robot_num].push_back(i_all_point_map[i_path_optimization[i_robot_num][i]]);

	current_robot_line[0] = path_optimization_real_xy[i_robot_num][0];
	current_robot_line[1] = path_optimization_real_xy[i_robot_num][1];

	for (crash_with_num = 1; crash_with_num < i_robot_num; crash_with_num++)
	{
		high_priority_robot_line[0] = path_optimization_real_xy[crash_with_num][0];
		high_priority_robot_line[1] = path_optimization_real_xy[crash_with_num][1];

		is_the_line_intersect[crash_with_num] = is_intersect(current_robot_line[0], current_robot_line[1], high_priority_robot_line[0], high_priority_robot_line[1]); //檢查是否相交

		if (is_the_line_intersect[crash_with_num])
		{
			have_intersection[crash_with_num] = find_intersection_point(current_robot_line[0], current_robot_line[1], high_priority_robot_line[0], high_priority_robot_line[1]); //解交點
			contact_distan_with_high_priority_robot[crash_with_num] = hypot(current_robot_line[0].x - have_intersection[crash_with_num].x, current_robot_line[0].y - have_intersection[crash_with_num].y);  //計算低優先權機器人與相交點的距離
			contact_distan_with_low_priority_robot[crash_with_num] = hypot(high_priority_robot_line[0].x - have_intersection[crash_with_num].x, high_priority_robot_line[0].y - have_intersection[crash_with_num].y); //計算高優先權機器人與相交點的距離

			if (
				(float)contact_distan_with_high_priority_robot[crash_with_num] / (float)contact_distan_with_low_priority_robot[crash_with_num] > 1.1
				|| contact_distan_with_high_priority_robot[crash_with_num] - contact_distan_with_low_priority_robot[crash_with_num]
				- (i_robot_num - crash_with_num) * 5 > 0 //再依據機器人的優先大小來決定是否要重新規劃路線 
				)
			{
				choose_path_num = 1;
				path_optimization[i_robot_num].clear();
				Path_Optimization_20181127(i_sca_image, i_all_point_map_original, i_robot_num, path_optimization[i_robot_num]);
			}
			else
			{
				choose_path_num = 0;
			}
		}
	}
//	choose_path_num = 0;
}

void CCollisionFreepathplanningDlg::binarization(Mat i_show_data, vector<vector<bool>>& o_sca_image2)
{
	vector<bool> sca_image1;

	unsigned char* i_show_data_data = i_show_data.data;

	for (int y = 0; y < i_show_data.rows - 1; y++)
	{
		for (int x = 0; x < i_show_data.cols - 1; x++)
		{
			int temp_of_image = i_show_data_data[y * i_show_data.rows + x];
			if (temp_of_image == 255)
				sca_image1.push_back(0);
			else
				sca_image1.push_back(1);
		}
		o_sca_image2.push_back(sca_image1);
		sca_image1.clear();
	}
/*
	for (int y = 0; y < i_show_data.rows - 1; y++)  //速度比較慢的取值方式
	{
		for (int x = 0; x < i_show_data.cols - 1; x++)
		{
			int temp_of_image = i_show_data.at<unsigned char>(y, x);
			if (temp_of_image == 255)
				sca_image1.push_back(0);
			else
				sca_image1.push_back(1);
		}
		o_sca_image2.push_back(sca_image1);
		sca_image1.clear();
	}*/

	for (int x = 0; x < i_show_data.cols - 1; x++)
	{
		sca_image1.push_back(0);
	}

	o_sca_image2.push_back(sca_image1);
	sca_image1.clear();

}

void CCollisionFreepathplanningDlg::find_coner(vector<vector<bool>> i_sca_image, vector <Point>& o_save_coner, int i_Interpolation)
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
					o_save_coner.push_back(Point(x, y));
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
				o_save_coner.push_back(Point(x, y));
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
					o_save_coner.push_back(Point(x, y));
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

void CCollisionFreepathplanningDlg::trans2Voronoi(vector<vector<bool>> i_sca_image, vector<Point> i_save_coner, double(&o_Data)[3000], int i_Interpolation2)
{
	unsigned int input_Data = 0;
	int jump_count = 0;
	jump_count = (i_sca_image.size() + 1) / i_Interpolation2;


	o_Data[0] = i_save_coner.size();
	for (input_Data = 1; input_Data < i_save_coner.size() + 1; input_Data++)
	{
		//		line(RGB_show_data, Point(save_coner[input_Data - 1].x, save_coner[input_Data - 1].y),Point(save_coner[input_Data - 1].x, save_coner[input_Data - 1].y), CV_RGB(0, 255, 0), 1);
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
				//				line(RGB_show_data, Point(i, j),Point(i, j), CV_RGB(200, 200, 0), 1);
			}
		}
	}
}

void CCollisionFreepathplanningDlg::boost_Voronoi_calculate(int x_boundary, int y_boundary, Point2d(&o_savepoint1)[1500], Point2d(&o_savepoint2)[1500], int& o_line_count)
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
void CCollisionFreepathplanningDlg::iterate_primary_edges1(const voronoi_diagram<double>& vd, Point2d(&o_savepoint1)[1500], Point2d(&o_savepoint2)[1500], int& o_line_count)
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

void CCollisionFreepathplanningDlg::Generalized_Voronoi(vector<vector<bool>> i_sca_image, Point2d i_savepoint1[1500], Point2d i_savepoint2[1500], int i_line_count, int& o_new_input_index, Point2d(&o_new_savepoint1)[1500], Point2d(&o_new_savepoint2)[1500])
{
	//	remove("廣義VD座標輸出.txt");
	//	fstream app_VD_output("廣義VD座標輸出.txt", ios::app);

	int cheak_point = 0, line_distant, new_input_index = 0;
	Point2d cutout_point[500] = { 0 };

	for (cheak_point = 0; cheak_point < i_line_count; cheak_point++)
	{
		if (i_savepoint1[cheak_point].x == 0 && i_savepoint1[cheak_point].y == 0 && i_savepoint2[cheak_point].x == 0 && i_savepoint2[cheak_point].y == 0)
			continue;

		if (i_savepoint1[cheak_point].x >= image_size - 1)  //拯救右邊界
			i_savepoint1[cheak_point].x = image_size - 1;
		if (i_savepoint2[cheak_point].x >= image_size - 1)
			i_savepoint2[cheak_point].x = image_size - 1;
		if (i_savepoint1[cheak_point].y >= image_size) //拯救下邊界
			i_savepoint1[cheak_point].y = image_size;
		if (i_savepoint2[cheak_point].y >= image_size)
			i_savepoint2[cheak_point].y = image_size;
		if (i_savepoint1[cheak_point].y <= 0)  //拯救上邊界
			i_savepoint1[cheak_point].y = 0;
		if (i_savepoint2[cheak_point].y <= 0)
			i_savepoint2[cheak_point].y = 0;

		line_distant = sqrt(pow(i_savepoint1[cheak_point].x - i_savepoint2[cheak_point].x, 2) + pow(i_savepoint1[cheak_point].y - i_savepoint2[cheak_point].y, 2));
		//		line_distant++;

		int cutout = 0;
		for (cutout = 0; cutout < line_distant + 1; cutout++)
		{
			//計算斜線上的點
			if (i_savepoint1[cheak_point].x >= i_savepoint2[cheak_point].x && i_savepoint1[cheak_point].y >= i_savepoint2[cheak_point].y)
			{
				cutout_point[cutout].x = i_savepoint2[cheak_point].x + (sqrt(pow((i_savepoint1[cheak_point].x - i_savepoint2[cheak_point].x), 2)) / (double)line_distant) * (double)cutout;
				cutout_point[cutout].y = i_savepoint2[cheak_point].y + (sqrt(pow((i_savepoint1[cheak_point].y - i_savepoint2[cheak_point].y), 2)) / (double)line_distant) * (double)cutout;
			}

			if (i_savepoint2[cheak_point].x >= i_savepoint1[cheak_point].x && i_savepoint1[cheak_point].y >= i_savepoint2[cheak_point].y)
			{
				cutout_point[cutout].x = i_savepoint1[cheak_point].x + (sqrt(pow((i_savepoint1[cheak_point].x - i_savepoint2[cheak_point].x), 2)) / (double)line_distant) * (double)cutout;
				cutout_point[cutout].y = i_savepoint1[cheak_point].y - (sqrt(pow((i_savepoint1[cheak_point].y - i_savepoint2[cheak_point].y), 2)) / (double)line_distant) * (double)cutout;
			}

			if (i_savepoint1[cheak_point].x >= i_savepoint2[cheak_point].x && i_savepoint2[cheak_point].y >= i_savepoint1[cheak_point].y)
			{
				cutout_point[cutout].x = i_savepoint2[cheak_point].x + (sqrt(pow((i_savepoint1[cheak_point].x - i_savepoint2[cheak_point].x), 2)) / (double)line_distant) * (double)cutout;
				cutout_point[cutout].y = i_savepoint2[cheak_point].y - (sqrt(pow((i_savepoint1[cheak_point].y - i_savepoint2[cheak_point].y), 2)) / (double)line_distant) * (double)cutout;
			}

			if (i_savepoint2[cheak_point].x >= i_savepoint1[cheak_point].x && i_savepoint2[cheak_point].y >= i_savepoint1[cheak_point].y)
			{
				cutout_point[cutout].x = i_savepoint1[cheak_point].x + (sqrt(pow((i_savepoint1[cheak_point].x - i_savepoint2[cheak_point].x), 2)) / (double)line_distant) * (double)cutout;
				cutout_point[cutout].y = i_savepoint1[cheak_point].y + (sqrt(pow((i_savepoint1[cheak_point].y - i_savepoint2[cheak_point].y), 2)) / (double)line_distant) * (double)cutout;
			}
		}

		if (cutout < 3)
		{
			if ((cutout_point[0].y > i_sca_image.size() - 1) || (cutout_point[0].x > i_sca_image[0].size() - 1) || (cutout_point[0].y < 1 || (cutout_point[0].x < 1)))
				continue;
			// 			if(_isnan(cutout_point[0].y))
			// 			{
			// 				continue;
			// 			}

			if (round(i_savepoint1[cheak_point].y > 87))
				i_savepoint1[cheak_point].y--;
			if (round(i_savepoint2[cheak_point].y > 87))
				i_savepoint2[cheak_point].y--;
			if (round(i_savepoint1[cheak_point].x > 87))
				i_savepoint1[cheak_point].x--;
			if (round(i_savepoint2[cheak_point].x > 87))
				i_savepoint2[cheak_point].x--;

			if (i_sca_image[round(i_savepoint1[cheak_point].y)][round(i_savepoint1[cheak_point].x)] == 0 &&/**/
				i_sca_image[round(i_savepoint2[cheak_point].y)][round(i_savepoint2[cheak_point].x)] == 0/**/)
			{
				o_new_savepoint1[new_input_index] = i_savepoint1[cheak_point];
				o_new_savepoint2[new_input_index] = i_savepoint2[cheak_point];
				new_input_index++;
			}
			continue;
		}

		for (int i = 0; i < cutout - 1; i++)
		{
			if ((cutout_point[i].y > i_sca_image.size() - 1) || (cutout_point[i].x > i_sca_image[0].size() - 1) || (cutout_point[i].y < 0 || (cutout_point[i].x < 0)))
				continue;

			if (i_sca_image[round(cutout_point[i].y)][round(cutout_point[i].x)] == 1
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
		memset(cutout_point, 0, sizeof(cutout_point));
	}
	o_new_input_index = new_input_index;  //輸出最後有多少廣義VD線段
}

void CCollisionFreepathplanningDlg::Match_point(int i_line_count, int i_new_input_index, Point2d(&io_new_savepoint1)[1500], Point2d(&io_new_savepoint2)[1500], float near_dis)
{
	int small_loop1, small_loop2, small_index = 0;
	double small_dis1, small_dis2, small_dis3;
	Point2d temp_point;
	vector <Point2d> center_point, center_point2;

	for (small_loop1 = 0; small_loop1 < i_line_count; small_loop1++)
	{
		for (small_loop2 = 0; small_loop2 < i_line_count; small_loop2++)
		{
			if (io_new_savepoint1[small_loop1].x == io_new_savepoint2[small_loop2].x &&
				io_new_savepoint1[small_loop1].y == io_new_savepoint2[small_loop2].y)
				continue;


			small_dis1 = sqrt(pow(io_new_savepoint1[small_loop1].x - io_new_savepoint2[small_loop2].x, 2) + pow(io_new_savepoint1[small_loop1].y - io_new_savepoint2[small_loop2].y, 2));
			if (small_dis1 < near_dis)
			{
				temp_point.x = (io_new_savepoint1[small_loop1].x + io_new_savepoint2[small_loop2].x) / 2;
				temp_point.y = (io_new_savepoint1[small_loop1].y + io_new_savepoint2[small_loop2].y) / 2;
				//io_new_savepoint2[small_loop2] = io_new_savepoint1[small_loop1];
				center_point.push_back(Point2d(temp_point.x, temp_point.y));
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

void CCollisionFreepathplanningDlg::Dijkstra_path_planning(int i_robot_num, Point2d i_robot_start[total_number], Point2d  i_robot_end[total_number], Point2d i_new_savepoint1[1500], Point2d i_new_savepoint2[1500], int i_new_input_index, vector <Point>& o_all_point_map, vector <Point2d>& o_all_point_map_original)
{
	vector <Point> CPoint_savepoint1;
	vector <Point> CPoint_savepoint2;
//	fstream app_path_num("path_num.txt", ios::app);
//	fstream app_path_num2("path_num2.txt", ios::app);


	for (int i = 0; i < i_new_input_index; i++)
	{
		CPoint_savepoint1.push_back(Point(round(i_new_savepoint1[i].x * 10), round(i_new_savepoint1[i].y * 10)));
		CPoint_savepoint2.push_back(Point(round(i_new_savepoint2[i].x * 10), round(i_new_savepoint2[i].y * 10)));
//		app_path_num2 << CPoint_savepoint1[i].x << ", " << CPoint_savepoint1[i].y << " 到 " << CPoint_savepoint2[i].x << ", " << CPoint_savepoint2[i].y << endl;
	}

	int loop1, loop2;
	bool onestime;
	int put_index = 0;

	//先丟入第一組，以免檢查時vector時沒東西
	o_all_point_map.push_back(CPoint_savepoint1[0]);
	o_all_point_map_original.push_back(i_new_savepoint1[0]);
//	app_path_num << o_all_point_map[put_index].x << ", " << o_all_point_map[put_index].y << " 第 " << put_index << endl;
	put_index++;

	//---------將CPoint_savepoint1與CPoint_savepoint2中所有不重複的點抓出來-----------------
	for (int i = 0; i < i_new_input_index; i++)
	{
		onestime = 1;
		for (unsigned int j = 0; j < o_all_point_map.size(); j++)
		{
			if (o_all_point_map[j] == CPoint_savepoint1[i])
			{
				onestime = 0;
				break;
			}
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
			{
				onestime = 0;
				break;
			}
		}
		if (onestime == 1)
		{
			o_all_point_map.push_back(CPoint_savepoint2[i]);
			o_all_point_map_original.push_back(i_new_savepoint2[i]); // 原始的也排列一次
//			app_path_num << o_all_point_map[put_index].x << ", " << o_all_point_map[put_index].y << " 第 " << put_index << endl;
			put_index++;
		}
	}
	//-----------------------------------------------------------------------------
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

	o_all_point_map.push_back(Point(i_robot_start[i_robot_num].x, i_robot_start[i_robot_num].y));
	CPoint_savepoint1.push_back(Point(i_robot_start[i_robot_num].x, i_robot_start[i_robot_num].y));
	o_all_point_map_original.push_back(Point2d(i_robot_start[i_robot_num].x / 10, i_robot_start[i_robot_num].y / 10));
	CPoint_savepoint2.push_back(Point(o_all_point_map[the_point_index[0]].x, o_all_point_map[the_point_index[0]].y));
	i_new_input_index++;
	put_index++;

	if (CPoint_savepoint1[i_new_input_index - 1] == CPoint_savepoint2[i_new_input_index - 1]) //為了解決機器人當下位置與VD點重合的狀況
	{
		CPoint_savepoint1[i_new_input_index - 1].x++;  //此刻機器人位置稍微位移
		o_all_point_map[o_all_point_map.size() - 1].x++;
	}


//	app_path_num2 << CPoint_savepoint1[i_new_input_index - 1].x << ", " << CPoint_savepoint1[i_new_input_index - 1].y << " 到 " << CPoint_savepoint2[i_new_input_index - 1].x << ", " << CPoint_savepoint2[i_new_input_index - 1].y << endl;

//		新的K-ShortestPath (YenTopKShortestPathsAlg) 20190710
	Graph my_graph(put_index, i_new_input_index, CPoint_savepoint1, CPoint_savepoint2, o_all_point_map);
	YenTopKShortestPathsAlg yenAlg(my_graph, my_graph.get_vertex(put_index - 1), my_graph.get_vertex(the_point_index[1]));

	int next_path = 0;
	while (yenAlg.has_next() && next_path < 9)
	{
		yenAlg.next()->MyOut(show_path[next_path], show_path_weight[next_path]);
		next_path++;
	}

// 	app_path_num.close();
// 	app_path_num2.close();
// 
// 	remove("path_num.txt");
// 	remove("path_num2.txt");
//    舊的dijkstra，於20190710拿除
// 	for (loop1 = 0; loop1 < i_new_input_index; loop1++)
// 	{
// 		for (loop2 = 0; loop2 < put_index; loop2++)
// 		{
// 
// 			if (CPoint_savepoint1[loop1] == o_all_point_map[loop2])
// 				savetemp_index1 = loop2;
// 
// 			if (CPoint_savepoint2[loop1] == o_all_point_map[loop2])
// 				savetemp_index2 = loop2;
// 
// 			if (w[loop1][loop2] == 0)
// 				w[loop1][loop2] = 50000;
// 
// 		}
// 
// 		w[savetemp_index1][savetemp_index2] = sqrt(pow(CPoint_savepoint1[loop1].x - CPoint_savepoint2[loop1].x, 2) + pow(CPoint_savepoint1[loop1].y - CPoint_savepoint2[loop1].y, 2));
// 		w[savetemp_index2][savetemp_index1] = w[savetemp_index1][savetemp_index2];
// 	}
// 
// 
// 
// 	dijkstra(put_index - 1, o_all_point_map.size());
// 	find_path(the_point_index[1]);
// 
// 
// 	if (check_path_change2.GetCheck() == TRUE)  //透過畫面上的check box來開關
// 	{
// 		int path_count = 0;
// 		int total_show_path = 0;
// 		int check_available_path_counter = 0;
// 		for (total_show_path = 1; total_show_path < show_path[0].size(); total_show_path++)
// 		{
// 			for (check_available_path_counter = 0; check_available_path_counter < o_all_point_map.size(); check_available_path_counter++)
// 			{
// 				if (w[show_path[0][total_show_path]][check_available_path_counter] < 50000)
// 					path_count++;
// 			}
// 
// 			if (path_count > 2)
// 				break;
// 
// 			path_count = 0;
// 		}
// 
// 		w[show_path[0][total_show_path]][show_path[0][total_show_path + 1]] = 50000;  //測試把最短路徑的第二節點到第三節點給設成無限大
// 		w[show_path[0][total_show_path + 1]][show_path[0][total_show_path]] = 50000;
// 		show_path[0].clear();
// 		dijkstra(put_index - 1, o_all_point_map.size());
// 		find_path(the_point_index[1]);
// 	}

}

void CCollisionFreepathplanningDlg::Control_Methods(int control_type, double i_rho, double i_alpha, double i_beta, double i_phi, double& o_vr, double& o_vl, int& o_state)
{
	//----------------------第二模式新增之判斷------------------------------------------------------------
	float P2_1 = 1.1737;
	float P2_2 = 1.4317;
	float P2_3 = 0.2422;
	float P2_4 = 0.4223;

	float Vt2 = i_rho * i_rho * P2_1 +
		i_alpha * i_alpha * P2_2 +
		i_alpha * i_phi * P2_4 +
		i_alpha * i_phi * P2_4 +
		i_phi * i_phi * P2_3;

	int rho_gain;
	int times = 0;
	float M1, M2, N1, N2;
	float alpha_gain;
	float beta_gain;
	//-----------------------------------------------------------------------------------------------------


	if (control_type == 1)
	{
		//原始線性控制
		o_vr = 2 * i_rho + 0.15 * (5 * i_alpha - 2 * (i_beta));
		o_vl = 2 * i_rho - 0.15 * (5 * i_alpha - 2 * (i_beta));
		//		o_vr = o_vr / 4;
		//		o_vl = o_vl / 4;


		//-------------------------normalize----------------------------------
		double temp = abs(o_vr) + abs(o_vl);
		o_vr = 100 * o_vr / temp;
		o_vl = 100 * o_vl / temp;
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
			o_vr = 200 * o_vr / temp;
			o_vl = 200 * o_vl / temp;
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

			o_vr = 0 * i_rho + 0.15 * (alpha_gain * i_alpha - 0 * (i_beta));
			o_vl = 0 * i_rho - 0.15 * (alpha_gain * i_alpha - 0 * (i_beta));
			o_vr = o_vr * 1.4;
			o_vl = o_vl * 1.4;

			break;

		case 2:  //直線Mode
				 // 			rho_gain = (4 * i_rho / 200) + (1 - (i_rho / 200));
				 // 			if (i_rho > 150) rho_gain = 4;
			rho_gain = 1.0331;
			o_vr = rho_gain * i_rho + 0.15 * (alpha_gain * i_alpha - 0 * (i_beta));
			o_vl = rho_gain * i_rho - 0.15 * (alpha_gain * i_alpha - 0 * (i_beta));
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
				N1 = (0.49 * CV_PI * sin(i_alpha) - sin(0.49 * CV_PI) * i_alpha) / (i_alpha * (0.49 * CV_PI - sin(0.49 * CV_PI)));

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

			alpha_gain = M1 * N1 * 1.2833 +
				M1 * N2 * 1.1022 +
				M2 * N1 * 1.2833 +
				M2 * N2 * 1.1022;

			beta_gain = -M1 * N1 * 0.0487 +
				-M1 * N2 * 0.0517 +
				-M2 * N1 * 0.0487 +
				-M2 * N2 * 0.0517;

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

void CCollisionFreepathplanningDlg::simulation_car(Mat& live_show, Point2d i_robot_start_point[total_number], double i_robot_zdir[total_number], int i_car_num, int i_carsize, int i_arrive_goal_times[total_number])
{
	//	fstream app_pos_m_sim("control_pos_m_sim.txt", ios::app);
	int total_arrive_times = 0;
	Point car[6];
	char num1[200];
	Point temp_path, detect_area[total_number][2] = { 0 };
	car_front_detector(i_robot_start_point, robot_zdir, detect_angle, i_carsize * (detecte_dis + 3) + 10, i_carsize, detect_area); //偵測前方是誰，並輸出繪圖用偵測範圍

	car[0] = Point(i_robot_start_point[i_car_num].x + i_carsize * cos(robot_zdir[i_car_num] + 0.7854), i_robot_start_point[i_car_num].y - i_carsize * sin(robot_zdir[i_car_num] + 0.7854));
	car[1] = Point(i_robot_start_point[i_car_num].x + i_carsize * cos(0.7854 - robot_zdir[i_car_num]), i_robot_start_point[i_car_num].y + i_carsize * sin(0.7854 - robot_zdir[i_car_num]));
	car[2] = Point(i_robot_start_point[i_car_num].x - i_carsize * cos(robot_zdir[i_car_num] + 0.7854), i_robot_start_point[i_car_num].y + i_carsize * sin(robot_zdir[i_car_num] + 0.7854));
	car[3] = Point(i_robot_start_point[i_car_num].x - i_carsize * cos(0.7854 - robot_zdir[i_car_num]), i_robot_start_point[i_car_num].y - i_carsize * sin(0.7854 - robot_zdir[i_car_num]));
	car[4] = Point(i_robot_start_point[i_car_num].x, i_robot_start_point[i_car_num].y);
	car[5] = Point(i_robot_start_point[i_car_num].x + i_carsize * 2 * cos(robot_zdir[i_car_num]), i_robot_start_point[i_car_num].y - i_carsize * 2 * sin(robot_zdir[i_car_num]));


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
	}
	else
	{
		circle(live_show, car[4], i_carsize - 2, CV_RGB(0, 255, 0), 2);
		// 		line(live_show, car[0], car[1], CV_RGB(0, 255, 0), 2);
		// 		line(live_show, car[1], car[2], CV_RGB(0, 255, 0), 2);
		// 		line(live_show, car[2], car[3], CV_RGB(0, 255, 0), 2);
		// 		line(live_show, car[3], car[0], CV_RGB(0, 255, 0), 2);
		line(live_show, car[4], car[5], CV_RGB(255, 150, 150), 2);
		if (detect_area[i_car_num][0].x != 0 && 1)  //顯示出偵測範圍以及半徑
		{
			line(live_show, car[4], detect_area[i_car_num][0], CV_RGB(255, 150, 150), 2);
			line(live_show, car[4], detect_area[i_car_num][1], CV_RGB(255, 150, 150), 2);
			double draw_x[2], draw_y[2], draw_angle[2];
			draw_x[0] = -(car[4].x - detect_area[i_car_num][0].x);
			draw_y[0] = -(car[4].y - detect_area[i_car_num][0].y);
			draw_x[1] = -(car[4].x - detect_area[i_car_num][1].x);
			draw_y[1] = -(car[4].y - detect_area[i_car_num][1].y);
			draw_angle[0] = atan2(draw_y[0], draw_x[0]);
			draw_angle[1] = atan2(draw_y[1], draw_x[1]);


			draw_angle[0] = draw_angle[0] * 180 / CV_PI;
			draw_angle[1] = draw_angle[1] * 180 / CV_PI;

			if (draw_x[1] < 0 && draw_y[1] < 0)
				draw_angle[1] = 360 + draw_angle[1];
			if (draw_x[1] >= 0 && draw_y[1] <= 0)
				draw_angle[1] = 360 + draw_angle[1];

			//line(live_show, detect_area[i_car_num][0], detect_area[i_car_num][1], CV_RGB(255, 150, 150), 1);

			ellipse(live_show, car[4], Size(i_carsize * (detecte_dis + 3) + 10, i_carsize * (detecte_dis + 3) + 10), 0, draw_angle[1], draw_angle[0], CV_RGB(255, 150, 150), 2, LINE_AA);
			//cvCircle(live_show, car[4], i_carsize * (detecte_dis + 2) + 10, CV_RGB(255, 150, 150), 2);
		}

		for (int i = all_robot_path[i_car_num].size(); i > 1; i--)
		{
			if (i == all_robot_path[i_car_num].size() - 500)
				break;
			if (i_car_num == 1)
				line(live_show, all_robot_path[i_car_num][i - 1], all_robot_path[i_car_num][i - 2], CV_RGB(255, 0, 0), 2, LINE_8);
			else if (i_car_num == 2 && total_number > 2)
				line(live_show, all_robot_path[i_car_num][i - 1], all_robot_path[i_car_num][i - 2], CV_RGB(100, 100, 100), 2, LINE_8);
			else if (i_car_num == 3 && total_number > 3)
				line(live_show, all_robot_path[i_car_num][i - 1], all_robot_path[i_car_num][i - 2], CV_RGB(220, 220, 0), 2, LINE_8);
			else if (i_car_num == 4 && total_number > 4)
				line(live_show, all_robot_path[i_car_num][i - 1], all_robot_path[i_car_num][i - 2], CV_RGB(170, 120, 50), 2, LINE_8);
			else if (i_car_num == 5 && total_number > 5)
				line(live_show, all_robot_path[i_car_num][i - 1], all_robot_path[i_car_num][i - 2], CV_RGB(180, 50, 100), 2, LINE_8);
			else if (i_car_num == 6 && total_number > 6)
				line(live_show, all_robot_path[i_car_num][i - 1], all_robot_path[i_car_num][i - 2], CV_RGB(50, 50, 120), 2, LINE_8);
			else if (i_car_num == 7 && total_number > 7)
				line(live_show, all_robot_path[i_car_num][i - 1], all_robot_path[i_car_num][i - 2], CV_RGB(160, 50, 170), 2, LINE_8);
			else if (i_car_num > 8 && total_number > 8)
				line(live_show, all_robot_path[i_car_num][i - 1], all_robot_path[i_car_num][i - 2], CV_RGB(230, 170, 100), 2, LINE_8);
		}
		sprintf_s(num1, "Robot %d", i_car_num);  //機器人編號

		putText(live_show, num1, Point(i_robot_start_point[i_car_num].x - 90, i_robot_start_point[i_car_num].y - i_carsize - 10), FONT_HERSHEY_PLAIN, 2, Scalar(70, 50, 50), 2, LINE_8);
		sprintf_s(num1, "%d", i_arrive_goal_times[i_car_num]); //個別到達次數
		putText(live_show, num1, Point(i_robot_start_point[i_car_num].x + 20, i_robot_start_point[i_car_num].y), FONT_HERSHEY_PLAIN, 1, Scalar(230, 20, 20), 1, LINE_8);

		sprintf_s(num1, "%.4f", show_path_efficiency[i_car_num]);//顯示路徑效率
		putText(live_show, num1, Point(i_robot_start_point[i_car_num].x + 20, i_robot_start_point[i_car_num].y + 30), FONT_HERSHEY_PLAIN, 1, Scalar(230, 20, 20), 1, LINE_8);

		sprintf_s(num1, "%d", run_times_counter);//顯示run_times_counter
		putText(live_show, num1, Point(800, 20), 1, 1, Scalar(230, 20, 20), 1, LINE_8);

	}



	// 	for (int i = 1; i < total_number; i++)
	// 	{
	// 		total_arrive_times += i_arrive_goal_times[i];//計算總共到達次數
	// 
	// 		sprintf_s(num1, "%d", car_front_num[i_car_num][i]);
	// 		putText(live_show, num1, Point(i_robot_start_point[i_car_num].x + i * 20, i_robot_start_point[i_car_num].y + 15), 1, 1, Scalar(230, 20, 20), 1, LINE_AA);
	// 
	// 
	// 	}

	// 	sprintf_s(num1, "total_arrive_times %d", total_arrive_times);
	// 	putText(live_show, num1, Point(650, 20), 1, 1, Scalar(230, 20, 20), 1, LINE_AA);
		//	app_pos_m_sim.close();
}

void CCollisionFreepathplanningDlg::car_crash_detector(Point2d i_robot_start_point[total_number], int i_carsize, int& io_crash_times, char io_crash_num[total_number][total_number])
{
	for (int i = 1; i < total_number; i++)
	{
		for (int j = i + 1; j < total_number; j++)
		{
			Point2d car_temp1 = i_robot_start_point[i], car_temp2 = i_robot_start_point[j];
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

void CCollisionFreepathplanningDlg::car_front_detector(Point2d i_robot_start_point[total_number], double i_robot_zdir[total_number], double i_detect_angle, int i_detect_radius, int i_carsize, Point(&o_detect_area)[total_number][2])
{
	for (int i = 1; i < total_number; i++)
	{
		o_detect_area[i][0] = Point(i_robot_start_point[i].x + i_detect_radius * cos(robot_zdir[i] + i_detect_angle), i_robot_start_point[i].y - i_detect_radius * sin(robot_zdir[i] + i_detect_angle));
		o_detect_area[i][1] = Point(i_robot_start_point[i].x + i_detect_radius * cos(i_detect_angle - robot_zdir[i]), i_robot_start_point[i].y + i_detect_radius * sin(i_detect_angle - robot_zdir[i]));

		Point2d car_temp1 = i_robot_start_point[i], car_direction_point = Point2d(i_robot_start_point[i].x + 20 * cos(robot_zdir[i]), i_robot_start_point[i].y - 20 * sin(robot_zdir[i]));
		double dx1 = car_direction_point.x - car_temp1.x;
		double dy1 = car_direction_point.y - car_temp1.y;

		for (int j = 1; j < total_number; j++)
		{
			if (i == j)
				continue;

			Point2d  car_temp2 = i_robot_start_point[j];
			double car_dis = sqrt(pow(car_temp1.x - car_temp2.x, 2) + pow(car_temp1.y - car_temp2.y, 2)) - i_carsize;
			double dx2 = car_temp2.x - car_temp1.x;
			double dy2 = car_temp2.y - car_temp1.y;
			double vector_angle = acos((dx1 * dx2 + dy1 * dy2) / sqrt((dx1 * dx1 + dy1 * dy1) * (dx2 * dx2 + dy2 * dy2) + 1e-10));

			if (car_dis < i_detect_radius + 1 && vector_angle < i_detect_angle + 0.01)
				car_front_num[i][j] = true;
			else
				car_front_num[i][j] = false;
		}
	}
}



void CCollisionFreepathplanningDlg::car_path_efficiency(Point2d i_robot_start_point, Point2d i_robot_end_point, int reentry_tag, int old_reentry_tag, int i_car_num)
{
	double ground_path_length = sqrt(pow(i_robot_start_point.x - i_robot_end_point.x, 2) + pow(i_robot_start_point.y - i_robot_end_point.y, 2));
	double real_path_length = 0;
	double temp_path = 0;
	double jump_caculate = 1;
	int i = 0;

	for (i = old_reentry_tag; i < reentry_tag - jump_caculate - 1; i += jump_caculate)
	{
		if (all_robot_path_double[i_car_num][i].x == all_robot_path_double[i_car_num][i + jump_caculate].x && all_robot_path_double[i_car_num][i].y == all_robot_path_double[i_car_num][i + jump_caculate].y)
			continue;

		temp_path = sqrt(pow(all_robot_path_double[i_car_num][i].x - all_robot_path_double[i_car_num][i + jump_caculate].x, 2) + pow(all_robot_path_double[i_car_num][i].y - all_robot_path_double[i_car_num][i + jump_caculate].y, 2)) + temp_path;
	}

	temp_path = sqrt(pow(all_robot_path_double[i_car_num][i].x - all_robot_path_double[i_car_num][reentry_tag - 1].x, 2) + pow(all_robot_path_double[i_car_num][i].y - all_robot_path_double[i_car_num][reentry_tag - 1].y, 2)) + temp_path + 4;

	show_path_efficiency[i_car_num] = /*ground_path_length /*/ temp_path;
	// 	if (show_path_efficiency[i_car_num] > 1)
// 		show_path_efficiency[i_car_num] = 1;
}


bool CCollisionFreepathplanningDlg::car_wait_mechanism(Point2d i_robot_start_point[total_number], double robot_zdir[total_number], int i_robot_num, vector <Point> i_all_point_map, vector <Point> i_path_optimization[total_number])
{
	//(paper 2)
	bool need_stop = false;
	int crash_with_num = 0, contact_distan_with_high_priority_robot[total_number] = { 0 }, contact_distan_with_low_priority_robot[total_number] = { 0 }, distan_of_current_robot[total_number] = { 0 }, distan_of_low_priority_robot[total_number] = { 0 }, distan_of_high_priority_robot[total_number] = { 0 };
	bool is_the_line_intersect[total_number] = { 0 };
	float intersection_angle[total_number] = { 0 };
	Point current_robot_line[2], low_priority_robot_line[2], high_priority_robot_line[2], have_intersection[total_number];

	for (int i = 1; i < total_number; i++)
	{
		if (i == i_robot_num)
			continue;

		if (car_front_num[i_robot_num][i] == true && i_robot_num > i)  //當低優先權偵測到高優先權的狀況
		{
			Point robot_line_translation[2];
			Point temp_point = current_robot_line[0] - high_priority_robot_line[0];

			need_stop = false;

			current_robot_line[0] = i_path_optimization[i_robot_num][0];
			current_robot_line[1] = i_path_optimization[i_robot_num][1];
			high_priority_robot_line[0] = i_path_optimization[i][0];
			high_priority_robot_line[1] = i_path_optimization[i][1];

			is_the_line_intersect[i_robot_num] = is_intersect(current_robot_line[0], current_robot_line[1], high_priority_robot_line[0], high_priority_robot_line[1]);

			if (is_the_line_intersect[i_robot_num])
			{
				have_intersection[i_robot_num] = find_intersection_point(current_robot_line[0], current_robot_line[1], high_priority_robot_line[0], high_priority_robot_line[1]); //解交點
				contact_distan_with_low_priority_robot[i_robot_num] = hypot(current_robot_line[0].x - have_intersection[i_robot_num].x, current_robot_line[0].y - have_intersection[i_robot_num].y);  //計算低優先權機器人與相交點的距離
				contact_distan_with_high_priority_robot[i] = hypot(high_priority_robot_line[0].x - have_intersection[i_robot_num].x, high_priority_robot_line[0].y - have_intersection[i_robot_num].y); //計算高優先權機器人與相交點的距離
				intersection_angle[i_robot_num] = find_intersection_angle(current_robot_line[0] - have_intersection[i_robot_num], high_priority_robot_line[0] - have_intersection[i_robot_num], (float)contact_distan_with_low_priority_robot[i_robot_num], (float)contact_distan_with_high_priority_robot[i]) * 180 / CV_PI;
				if (intersection_angle[i_robot_num] < 120)
					need_stop = true;
			}
		}


		if (car_front_num[i_robot_num][i] == true && i_robot_num < i)  //當高優先權偵測到低優先權的狀況
		{
			Point robot_line_translation[2];
			Point temp_point = current_robot_line[0] - low_priority_robot_line[0];

			current_robot_line[0] = i_path_optimization[i_robot_num][0];
			current_robot_line[1] = i_path_optimization[i_robot_num][1];
			low_priority_robot_line[0] = i_path_optimization[i][0];
			low_priority_robot_line[1] = i_path_optimization[i][1];
			robot_line_translation[0] = low_priority_robot_line[0] + temp_point;
			robot_line_translation[1] = low_priority_robot_line[1] + temp_point;

			is_the_line_intersect[i] = is_intersect(current_robot_line[0], current_robot_line[1], low_priority_robot_line[0], low_priority_robot_line[1]);

			if (is_the_line_intersect[i])
			{
				have_intersection[i] = find_intersection_point(current_robot_line[0], current_robot_line[1], low_priority_robot_line[0], low_priority_robot_line[1]); //解交點
				contact_distan_with_high_priority_robot[i_robot_num] = hypot(current_robot_line[0].x - have_intersection[i].x, current_robot_line[0].y - have_intersection[i].y);  //計算高優先權機器人與相交點的距離
				contact_distan_with_low_priority_robot[i] = hypot(low_priority_robot_line[0].x - have_intersection[i].x, low_priority_robot_line[0].y - have_intersection[i].y); //計算低優先權機器人與相交點的距離
				intersection_angle[i] = find_intersection_angle(current_robot_line[0] - have_intersection[i], low_priority_robot_line[0] - have_intersection[i], (float)contact_distan_with_high_priority_robot[i_robot_num], (float)contact_distan_with_low_priority_robot[i]) * 180 / CV_PI;
				if (intersection_angle[i] > 120)
					need_stop = true;
			}
			else if (car_front_num[i][i_robot_num] == true)
			{
				distan_of_current_robot[i_robot_num] = hypot(current_robot_line[0].x - current_robot_line[1].x, current_robot_line[0].y - current_robot_line[1].y);  //計算高優先權機器人下一刻路徑規劃的長度
				distan_of_low_priority_robot[i] = hypot(robot_line_translation[0].x - robot_line_translation[1].x, robot_line_translation[0].y - robot_line_translation[1].y); //計算低優先權機器人下一刻路徑規劃的長度
				intersection_angle[i] = find_intersection_angle(current_robot_line[1] - current_robot_line[0], robot_line_translation[1] - robot_line_translation[0], (float)distan_of_current_robot[i_robot_num], (float)distan_of_low_priority_robot[i]) * 180 / CV_PI;
				if (intersection_angle[i] > 120)
					need_stop = true;
			}

			if (car_front_num[i][i_robot_num] == false)
				need_stop = true;
		}
	}

	return need_stop;
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
	// 	Sleep(200);
	//	CDialogEx::OnOK();
}


