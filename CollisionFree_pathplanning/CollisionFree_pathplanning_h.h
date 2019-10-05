

/* this ALWAYS GENERATED file contains the definitions for the interfaces */


 /* File created by MIDL compiler version 8.01.0622 */
/* at Tue Jan 19 11:14:07 2038
 */
/* Compiler settings for CollisionFreepathplanning.idl:
    Oicf, W1, Zp8, env=Win64 (32b run), target_arch=AMD64 8.01.0622 
    protocol : all , ms_ext, c_ext, robust
    error checks: allocation ref bounds_check enum stub_data 
    VC __declspec() decoration level: 
         __declspec(uuid()), __declspec(selectany), __declspec(novtable)
         DECLSPEC_UUID(), MIDL_INTERFACE()
*/
/* @@MIDL_FILE_HEADING(  ) */



/* verify that the <rpcndr.h> version is high enough to compile this file*/
#ifndef __REQUIRED_RPCNDR_H_VERSION__
#define __REQUIRED_RPCNDR_H_VERSION__ 500
#endif

#include "rpc.h"
#include "rpcndr.h"

#ifndef __RPCNDR_H_VERSION__
#error this stub requires an updated version of <rpcndr.h>
#endif /* __RPCNDR_H_VERSION__ */


#ifndef __CollisionFree_pathplanning_h_h__
#define __CollisionFree_pathplanning_h_h__

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

/* Forward Declarations */ 

#ifndef __ICollisionFree_pathplanning_FWD_DEFINED__
#define __ICollisionFree_pathplanning_FWD_DEFINED__
typedef interface ICollisionFree_pathplanning ICollisionFree_pathplanning;

#endif 	/* __ICollisionFree_pathplanning_FWD_DEFINED__ */


#ifndef __CollisionFree_pathplanning_FWD_DEFINED__
#define __CollisionFree_pathplanning_FWD_DEFINED__

#ifdef __cplusplus
typedef class CollisionFree_pathplanning CollisionFree_pathplanning;
#else
typedef struct CollisionFree_pathplanning CollisionFree_pathplanning;
#endif /* __cplusplus */

#endif 	/* __CollisionFree_pathplanning_FWD_DEFINED__ */


#ifdef __cplusplus
extern "C"{
#endif 



#ifndef __CollisionFree_pathplanning_LIBRARY_DEFINED__
#define __CollisionFree_pathplanning_LIBRARY_DEFINED__

/* library CollisionFree_pathplanning */
/* [version][uuid] */ 


EXTERN_C const IID LIBID_CollisionFree_pathplanning;

#ifndef __ICollisionFree_pathplanning_DISPINTERFACE_DEFINED__
#define __ICollisionFree_pathplanning_DISPINTERFACE_DEFINED__

/* dispinterface ICollisionFree_pathplanning */
/* [uuid] */ 


EXTERN_C const IID DIID_ICollisionFree_pathplanning;

#if defined(__cplusplus) && !defined(CINTERFACE)

    MIDL_INTERFACE("7d1d0d09-f93c-4c2f-8b7a-4d0b86428d74")
    ICollisionFree_pathplanning : public IDispatch
    {
    };
    
#else 	/* C style interface */

    typedef struct ICollisionFree_pathplanningVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE *QueryInterface )( 
            ICollisionFree_pathplanning * This,
            /* [in] */ REFIID riid,
            /* [annotation][iid_is][out] */ 
            _COM_Outptr_  void **ppvObject);
        
        ULONG ( STDMETHODCALLTYPE *AddRef )( 
            ICollisionFree_pathplanning * This);
        
        ULONG ( STDMETHODCALLTYPE *Release )( 
            ICollisionFree_pathplanning * This);
        
        HRESULT ( STDMETHODCALLTYPE *GetTypeInfoCount )( 
            ICollisionFree_pathplanning * This,
            /* [out] */ UINT *pctinfo);
        
        HRESULT ( STDMETHODCALLTYPE *GetTypeInfo )( 
            ICollisionFree_pathplanning * This,
            /* [in] */ UINT iTInfo,
            /* [in] */ LCID lcid,
            /* [out] */ ITypeInfo **ppTInfo);
        
        HRESULT ( STDMETHODCALLTYPE *GetIDsOfNames )( 
            ICollisionFree_pathplanning * This,
            /* [in] */ REFIID riid,
            /* [size_is][in] */ LPOLESTR *rgszNames,
            /* [range][in] */ UINT cNames,
            /* [in] */ LCID lcid,
            /* [size_is][out] */ DISPID *rgDispId);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE *Invoke )( 
            ICollisionFree_pathplanning * This,
            /* [annotation][in] */ 
            _In_  DISPID dispIdMember,
            /* [annotation][in] */ 
            _In_  REFIID riid,
            /* [annotation][in] */ 
            _In_  LCID lcid,
            /* [annotation][in] */ 
            _In_  WORD wFlags,
            /* [annotation][out][in] */ 
            _In_  DISPPARAMS *pDispParams,
            /* [annotation][out] */ 
            _Out_opt_  VARIANT *pVarResult,
            /* [annotation][out] */ 
            _Out_opt_  EXCEPINFO *pExcepInfo,
            /* [annotation][out] */ 
            _Out_opt_  UINT *puArgErr);
        
        END_INTERFACE
    } ICollisionFree_pathplanningVtbl;

    interface ICollisionFree_pathplanning
    {
        CONST_VTBL struct ICollisionFree_pathplanningVtbl *lpVtbl;
    };

    

#ifdef COBJMACROS


#define ICollisionFree_pathplanning_QueryInterface(This,riid,ppvObject)	\
    ( (This)->lpVtbl -> QueryInterface(This,riid,ppvObject) ) 

#define ICollisionFree_pathplanning_AddRef(This)	\
    ( (This)->lpVtbl -> AddRef(This) ) 

#define ICollisionFree_pathplanning_Release(This)	\
    ( (This)->lpVtbl -> Release(This) ) 


#define ICollisionFree_pathplanning_GetTypeInfoCount(This,pctinfo)	\
    ( (This)->lpVtbl -> GetTypeInfoCount(This,pctinfo) ) 

#define ICollisionFree_pathplanning_GetTypeInfo(This,iTInfo,lcid,ppTInfo)	\
    ( (This)->lpVtbl -> GetTypeInfo(This,iTInfo,lcid,ppTInfo) ) 

#define ICollisionFree_pathplanning_GetIDsOfNames(This,riid,rgszNames,cNames,lcid,rgDispId)	\
    ( (This)->lpVtbl -> GetIDsOfNames(This,riid,rgszNames,cNames,lcid,rgDispId) ) 

#define ICollisionFree_pathplanning_Invoke(This,dispIdMember,riid,lcid,wFlags,pDispParams,pVarResult,pExcepInfo,puArgErr)	\
    ( (This)->lpVtbl -> Invoke(This,dispIdMember,riid,lcid,wFlags,pDispParams,pVarResult,pExcepInfo,puArgErr) ) 

#endif /* COBJMACROS */


#endif 	/* C style interface */


#endif 	/* __ICollisionFree_pathplanning_DISPINTERFACE_DEFINED__ */


EXTERN_C const CLSID CLSID_CollisionFree_pathplanning;

#ifdef __cplusplus

class DECLSPEC_UUID("36d70757-6cd3-4745-9cc5-285035f97420")
CollisionFree_pathplanning;
#endif
#endif /* __CollisionFree_pathplanning_LIBRARY_DEFINED__ */

/* Additional Prototypes for ALL interfaces */

/* end of Additional Prototypes */

#ifdef __cplusplus
}
#endif

#endif


