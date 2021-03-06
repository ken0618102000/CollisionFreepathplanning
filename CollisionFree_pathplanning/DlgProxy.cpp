
// DlgProxy.cpp: 實作檔案
//

#include "stdafx.h"
#include "CollisionFree_pathplanning.h"
#include "DlgProxy.h"
#include "CollisionFree_pathplanningDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CCollisionFreepathplanningDlgAutoProxy

IMPLEMENT_DYNCREATE(CCollisionFreepathplanningDlgAutoProxy, CCmdTarget)

CCollisionFreepathplanningDlgAutoProxy::CCollisionFreepathplanningDlgAutoProxy()
{
	EnableAutomation();

	// 若要保持應用程式執行的時間與自動化一樣長
	//	建構函式必須呼叫 AfxOleLockApp。
	AfxOleLockApp();

	// 經由應用程式主視窗指標，取得對對話方塊的存取。
	// 將 Proxy 的內部指標設定為指向對話方塊，
	// 並設定指向此 Proxy 的
	// 對話方塊返回指標。
	ASSERT_VALID(AfxGetApp()->m_pMainWnd);
	if (AfxGetApp()->m_pMainWnd)
	{
		ASSERT_KINDOF(CCollisionFreepathplanningDlg, AfxGetApp()->m_pMainWnd);
		if (AfxGetApp()->m_pMainWnd->IsKindOf(RUNTIME_CLASS(CCollisionFreepathplanningDlg)))
		{
			m_pDialog = reinterpret_cast<CCollisionFreepathplanningDlg*>(AfxGetApp()->m_pMainWnd);
			m_pDialog->m_pAutoProxy = this;
		}
	}
}

CCollisionFreepathplanningDlgAutoProxy::~CCollisionFreepathplanningDlgAutoProxy()
{
	// 若要在使用 Automation 建立了所有物件之後結束應用程式，
	//	 解構函式必須呼叫 AfxOleUnlockApp。
	// 此外，這會摧毀主對話方塊
	if (m_pDialog != nullptr)
		m_pDialog->m_pAutoProxy = nullptr;
	AfxOleUnlockApp();
}

void CCollisionFreepathplanningDlgAutoProxy::OnFinalRelease()
{
	// 當釋放 Automation 物件最後的參考時，
	// 會呼叫 OnFinalRelease。基底類別會自動
	// 刪除物件。呼叫基底類別前，請先加入您物件所需的額外清除 (Cleanup)
	// 程式碼。

	CCmdTarget::OnFinalRelease();
}

BEGIN_MESSAGE_MAP(CCollisionFreepathplanningDlgAutoProxy, CCmdTarget)
END_MESSAGE_MAP()

BEGIN_DISPATCH_MAP(CCollisionFreepathplanningDlgAutoProxy, CCmdTarget)
END_DISPATCH_MAP()

// 備註: 我們新增了對 IID_ICollisionFree_pathplanning 的支援以支援型別安全繫結
//  以便從 VBA 支援類型安全繫結。此 IID 必須與 .IDL 檔中，
// 附加至分配介面 (Dispinterface) 的 GUID 相符。

// {7d1d0d09-f93c-4c2f-8b7a-4d0b86428d74}
static const IID IID_ICollisionFree_pathplanning =
{0x7d1d0d09,0xf93c,0x4c2f,{0x8b,0x7a,0x4d,0x0b,0x86,0x42,0x8d,0x74}};

BEGIN_INTERFACE_MAP(CCollisionFreepathplanningDlgAutoProxy, CCmdTarget)
	INTERFACE_PART(CCollisionFreepathplanningDlgAutoProxy, IID_ICollisionFree_pathplanning, Dispatch)
END_INTERFACE_MAP()

// 在此專案的 StdAfx.h 中定義 IMPLEMENT_OLECREATE2 巨集
// {36d70757-6cd3-4745-9cc5-285035f97420}
IMPLEMENT_OLECREATE2(CCollisionFreepathplanningDlgAutoProxy, "CollisionFree_pathplanning.Application", 0x36d70757,0x6cd3,0x4745,0x9c,0xc5,0x28,0x50,0x35,0xf9,0x74,0x20)


// CCollisionFreepathplanningDlgAutoProxy 訊息處理常式
