
// DlgProxy.h: 標頭檔
//

#pragma once

class CCollisionFreepathplanningDlg;


// CCollisionFreepathplanningDlgAutoProxy 命令目標

class CCollisionFreepathplanningDlgAutoProxy : public CCmdTarget
{
	DECLARE_DYNCREATE(CCollisionFreepathplanningDlgAutoProxy)

	CCollisionFreepathplanningDlgAutoProxy();           // 動態建立所使用的保護建構函式

// 屬性
public:
	CCollisionFreepathplanningDlg* m_pDialog;

// 作業
public:

// 覆寫
	public:
	virtual void OnFinalRelease();

// 程式碼實作
protected:
	virtual ~CCollisionFreepathplanningDlgAutoProxy();

	// 產生的訊息對應函式

	DECLARE_MESSAGE_MAP()
	DECLARE_OLECREATE(CCollisionFreepathplanningDlgAutoProxy)

	// 產生的 OLE 分派對應函式

	DECLARE_DISPATCH_MAP()
	DECLARE_INTERFACE_MAP()
};

