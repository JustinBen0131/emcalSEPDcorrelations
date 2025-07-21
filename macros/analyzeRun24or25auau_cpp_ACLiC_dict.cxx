// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIUsersdIpatsfan753dIDesktopdIauauAnalysisdIemcalSEPDcorrelationsdImacrosdIanalyzeRun24or25auau_cpp_ACLiC_dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp"

// Header files passed via #pragma extra_include

namespace term {
   namespace ROOTDict {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *term_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("term", 0 /*version*/, "analyzeRun24or25auau.cpp", 102,
                     ::ROOT::Internal::DefineBehavior((void*)nullptr,(void*)nullptr),
                     &term_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *term_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace log {
   namespace ROOTDict {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *log_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("log", 0 /*version*/, "analyzeRun24or25auau.cpp", 111,
                     ::ROOT::Internal::DefineBehavior((void*)nullptr,(void*)nullptr),
                     &log_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *log_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace ROOT {
   static TClass *CutKey_Dictionary();
   static void CutKey_TClassManip(TClass*);
   static void *new_CutKey(void *p = nullptr);
   static void *newArray_CutKey(Long_t size, void *p);
   static void delete_CutKey(void *p);
   static void deleteArray_CutKey(void *p);
   static void destruct_CutKey(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CutKey*)
   {
      ::CutKey *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::CutKey));
      static ::ROOT::TGenericClassInfo 
         instance("CutKey", "analyzeRun24or25auau.cpp", 205,
                  typeid(::CutKey), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &CutKey_Dictionary, isa_proxy, 4,
                  sizeof(::CutKey) );
      instance.SetNew(&new_CutKey);
      instance.SetNewArray(&newArray_CutKey);
      instance.SetDelete(&delete_CutKey);
      instance.SetDeleteArray(&deleteArray_CutKey);
      instance.SetDestructor(&destruct_CutKey);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CutKey*)
   {
      return GenerateInitInstanceLocal(static_cast<::CutKey*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::CutKey*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CutKey_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::CutKey*>(nullptr))->GetClass();
      CutKey_TClassManip(theClass);
   return theClass;
   }

   static void CutKey_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *NSCachelEMapPairgR_Dictionary();
   static void NSCachelEMapPairgR_TClassManip(TClass*);
   static void *new_NSCachelEMapPairgR(void *p = nullptr);
   static void *newArray_NSCachelEMapPairgR(Long_t size, void *p);
   static void delete_NSCachelEMapPairgR(void *p);
   static void deleteArray_NSCachelEMapPairgR(void *p);
   static void destruct_NSCachelEMapPairgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NSCache<MapPair>*)
   {
      ::NSCache<MapPair> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::NSCache<MapPair>));
      static ::ROOT::TGenericClassInfo 
         instance("NSCache<MapPair>", "analyzeRun24or25auau.cpp", 246,
                  typeid(::NSCache<MapPair>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &NSCachelEMapPairgR_Dictionary, isa_proxy, 4,
                  sizeof(::NSCache<MapPair>) );
      instance.SetNew(&new_NSCachelEMapPairgR);
      instance.SetNewArray(&newArray_NSCachelEMapPairgR);
      instance.SetDelete(&delete_NSCachelEMapPairgR);
      instance.SetDeleteArray(&deleteArray_NSCachelEMapPairgR);
      instance.SetDestructor(&destruct_NSCachelEMapPairgR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NSCache<MapPair>*)
   {
      return GenerateInitInstanceLocal(static_cast<::NSCache<MapPair>*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::NSCache<MapPair>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *NSCachelEMapPairgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::NSCache<MapPair>*>(nullptr))->GetClass();
      NSCachelEMapPairgR_TClassManip(theClass);
   return theClass;
   }

   static void NSCachelEMapPairgR_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *NSCachelEMapPairgRcLcLH_Dictionary();
   static void NSCachelEMapPairgRcLcLH_TClassManip(TClass*);
   static void *new_NSCachelEMapPairgRcLcLH(void *p = nullptr);
   static void *newArray_NSCachelEMapPairgRcLcLH(Long_t size, void *p);
   static void delete_NSCachelEMapPairgRcLcLH(void *p);
   static void deleteArray_NSCachelEMapPairgRcLcLH(void *p);
   static void destruct_NSCachelEMapPairgRcLcLH(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NSCache<MapPair>::H*)
   {
      ::NSCache<MapPair>::H *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::NSCache<MapPair>::H));
      static ::ROOT::TGenericClassInfo 
         instance("NSCache<MapPair>::H", "analyzeRun24or25auau.cpp", 251,
                  typeid(::NSCache<MapPair>::H), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &NSCachelEMapPairgRcLcLH_Dictionary, isa_proxy, 4,
                  sizeof(::NSCache<MapPair>::H) );
      instance.SetNew(&new_NSCachelEMapPairgRcLcLH);
      instance.SetNewArray(&newArray_NSCachelEMapPairgRcLcLH);
      instance.SetDelete(&delete_NSCachelEMapPairgRcLcLH);
      instance.SetDeleteArray(&deleteArray_NSCachelEMapPairgRcLcLH);
      instance.SetDestructor(&destruct_NSCachelEMapPairgRcLcLH);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NSCache<MapPair>::H*)
   {
      return GenerateInitInstanceLocal(static_cast<::NSCache<MapPair>::H*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::NSCache<MapPair>::H*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *NSCachelEMapPairgRcLcLH_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::NSCache<MapPair>::H*>(nullptr))->GetClass();
      NSCachelEMapPairgRcLcLH_TClassManip(theClass);
   return theClass;
   }

   static void NSCachelEMapPairgRcLcLH_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *MapPair_Dictionary();
   static void MapPair_TClassManip(TClass*);
   static void *new_MapPair(void *p = nullptr);
   static void *newArray_MapPair(Long_t size, void *p);
   static void delete_MapPair(void *p);
   static void deleteArray_MapPair(void *p);
   static void destruct_MapPair(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MapPair*)
   {
      ::MapPair *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::MapPair));
      static ::ROOT::TGenericClassInfo 
         instance("MapPair", "analyzeRun24or25auau.cpp", 259,
                  typeid(::MapPair), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &MapPair_Dictionary, isa_proxy, 4,
                  sizeof(::MapPair) );
      instance.SetNew(&new_MapPair);
      instance.SetNewArray(&newArray_MapPair);
      instance.SetDelete(&delete_MapPair);
      instance.SetDeleteArray(&deleteArray_MapPair);
      instance.SetDestructor(&destruct_MapPair);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MapPair*)
   {
      return GenerateInitInstanceLocal(static_cast<::MapPair*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::MapPair*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *MapPair_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::MapPair*>(nullptr))->GetClass();
      MapPair_TClassManip(theClass);
   return theClass;
   }

   static void MapPair_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *QA_Dictionary();
   static void QA_TClassManip(TClass*);
   static void delete_QA(void *p);
   static void deleteArray_QA(void *p);
   static void destruct_QA(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::QA*)
   {
      ::QA *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::QA));
      static ::ROOT::TGenericClassInfo 
         instance("QA", "analyzeRun24or25auau.cpp", 266,
                  typeid(::QA), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &QA_Dictionary, isa_proxy, 4,
                  sizeof(::QA) );
      instance.SetDelete(&delete_QA);
      instance.SetDeleteArray(&deleteArray_QA);
      instance.SetDestructor(&destruct_QA);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::QA*)
   {
      return GenerateInitInstanceLocal(static_cast<::QA*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::QA*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *QA_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::QA*>(nullptr))->GetClass();
      QA_TClassManip(theClass);
   return theClass;
   }

   static void QA_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pi0QA_Dictionary();
   static void Pi0QA_TClassManip(TClass*);
   static void delete_Pi0QA(void *p);
   static void deleteArray_Pi0QA(void *p);
   static void destruct_Pi0QA(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pi0QA*)
   {
      ::Pi0QA *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pi0QA));
      static ::ROOT::TGenericClassInfo 
         instance("Pi0QA", "analyzeRun24or25auau.cpp", 286,
                  typeid(::Pi0QA), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pi0QA_Dictionary, isa_proxy, 4,
                  sizeof(::Pi0QA) );
      instance.SetDelete(&delete_Pi0QA);
      instance.SetDeleteArray(&deleteArray_Pi0QA);
      instance.SetDestructor(&destruct_Pi0QA);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pi0QA*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pi0QA*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pi0QA*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pi0QA_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pi0QA*>(nullptr))->GetClass();
      Pi0QA_TClassManip(theClass);
   return theClass;
   }

   static void Pi0QA_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pi0QAcLcLFitInfo_Dictionary();
   static void Pi0QAcLcLFitInfo_TClassManip(TClass*);
   static void *new_Pi0QAcLcLFitInfo(void *p = nullptr);
   static void *newArray_Pi0QAcLcLFitInfo(Long_t size, void *p);
   static void delete_Pi0QAcLcLFitInfo(void *p);
   static void deleteArray_Pi0QAcLcLFitInfo(void *p);
   static void destruct_Pi0QAcLcLFitInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pi0QA::FitInfo*)
   {
      ::Pi0QA::FitInfo *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pi0QA::FitInfo));
      static ::ROOT::TGenericClassInfo 
         instance("Pi0QA::FitInfo", "analyzeRun24or25auau.cpp", 293,
                  typeid(::Pi0QA::FitInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pi0QAcLcLFitInfo_Dictionary, isa_proxy, 4,
                  sizeof(::Pi0QA::FitInfo) );
      instance.SetNew(&new_Pi0QAcLcLFitInfo);
      instance.SetNewArray(&newArray_Pi0QAcLcLFitInfo);
      instance.SetDelete(&delete_Pi0QAcLcLFitInfo);
      instance.SetDeleteArray(&deleteArray_Pi0QAcLcLFitInfo);
      instance.SetDestructor(&destruct_Pi0QAcLcLFitInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pi0QA::FitInfo*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pi0QA::FitInfo*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pi0QA::FitInfo*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pi0QAcLcLFitInfo_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pi0QA::FitInfo*>(nullptr))->GetClass();
      Pi0QAcLcLFitInfo_TClassManip(theClass);
   return theClass;
   }

   static void Pi0QAcLcLFitInfo_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pi0QAcLcLFitPair_Dictionary();
   static void Pi0QAcLcLFitPair_TClassManip(TClass*);
   static void *new_Pi0QAcLcLFitPair(void *p = nullptr);
   static void *newArray_Pi0QAcLcLFitPair(Long_t size, void *p);
   static void delete_Pi0QAcLcLFitPair(void *p);
   static void deleteArray_Pi0QAcLcLFitPair(void *p);
   static void destruct_Pi0QAcLcLFitPair(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pi0QA::FitPair*)
   {
      ::Pi0QA::FitPair *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pi0QA::FitPair));
      static ::ROOT::TGenericClassInfo 
         instance("Pi0QA::FitPair", "analyzeRun24or25auau.cpp", 299,
                  typeid(::Pi0QA::FitPair), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pi0QAcLcLFitPair_Dictionary, isa_proxy, 4,
                  sizeof(::Pi0QA::FitPair) );
      instance.SetNew(&new_Pi0QAcLcLFitPair);
      instance.SetNewArray(&newArray_Pi0QAcLcLFitPair);
      instance.SetDelete(&delete_Pi0QAcLcLFitPair);
      instance.SetDeleteArray(&deleteArray_Pi0QAcLcLFitPair);
      instance.SetDestructor(&destruct_Pi0QAcLcLFitPair);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pi0QA::FitPair*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pi0QA::FitPair*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pi0QA::FitPair*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pi0QAcLcLFitPair_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pi0QA::FitPair*>(nullptr))->GetClass();
      Pi0QAcLcLFitPair_TClassManip(theClass);
   return theClass;
   }

   static void Pi0QAcLcLFitPair_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *NSPair_Dictionary();
   static void NSPair_TClassManip(TClass*);
   static void *new_NSPair(void *p = nullptr);
   static void *newArray_NSPair(Long_t size, void *p);
   static void delete_NSPair(void *p);
   static void deleteArray_NSPair(void *p);
   static void destruct_NSPair(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NSPair*)
   {
      ::NSPair *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::NSPair));
      static ::ROOT::TGenericClassInfo 
         instance("NSPair", "analyzeRun24or25auau.cpp", 1274,
                  typeid(::NSPair), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &NSPair_Dictionary, isa_proxy, 4,
                  sizeof(::NSPair) );
      instance.SetNew(&new_NSPair);
      instance.SetNewArray(&newArray_NSPair);
      instance.SetDelete(&delete_NSPair);
      instance.SetDeleteArray(&deleteArray_NSPair);
      instance.SetDestructor(&destruct_NSPair);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NSPair*)
   {
      return GenerateInitInstanceLocal(static_cast<::NSPair*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::NSPair*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *NSPair_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::NSPair*>(nullptr))->GetClass();
      NSPair_TClassManip(theClass);
   return theClass;
   }

   static void NSPair_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *CorrQA_Dictionary();
   static void CorrQA_TClassManip(TClass*);
   static void delete_CorrQA(void *p);
   static void deleteArray_CorrQA(void *p);
   static void destruct_CorrQA(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CorrQA*)
   {
      ::CorrQA *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::CorrQA));
      static ::ROOT::TGenericClassInfo 
         instance("CorrQA", "analyzeRun24or25auau.cpp", 1277,
                  typeid(::CorrQA), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &CorrQA_Dictionary, isa_proxy, 4,
                  sizeof(::CorrQA) );
      instance.SetDelete(&delete_CorrQA);
      instance.SetDeleteArray(&deleteArray_CorrQA);
      instance.SetDestructor(&destruct_CorrQA);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CorrQA*)
   {
      return GenerateInitInstanceLocal(static_cast<::CorrQA*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::CorrQA*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CorrQA_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::CorrQA*>(nullptr))->GetClass();
      CorrQA_TClassManip(theClass);
   return theClass;
   }

   static void CorrQA_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *EmcalQA_Dictionary();
   static void EmcalQA_TClassManip(TClass*);
   static void delete_EmcalQA(void *p);
   static void deleteArray_EmcalQA(void *p);
   static void destruct_EmcalQA(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EmcalQA*)
   {
      ::EmcalQA *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::EmcalQA));
      static ::ROOT::TGenericClassInfo 
         instance("EmcalQA", "analyzeRun24or25auau.cpp", 1656,
                  typeid(::EmcalQA), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &EmcalQA_Dictionary, isa_proxy, 4,
                  sizeof(::EmcalQA) );
      instance.SetDelete(&delete_EmcalQA);
      instance.SetDeleteArray(&deleteArray_EmcalQA);
      instance.SetDestructor(&destruct_EmcalQA);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EmcalQA*)
   {
      return GenerateInitInstanceLocal(static_cast<::EmcalQA*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::EmcalQA*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *EmcalQA_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::EmcalQA*>(nullptr))->GetClass();
      EmcalQA_TClassManip(theClass);
   return theClass;
   }

   static void EmcalQA_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *HcalQA_Dictionary();
   static void HcalQA_TClassManip(TClass*);
   static void delete_HcalQA(void *p);
   static void deleteArray_HcalQA(void *p);
   static void destruct_HcalQA(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HcalQA*)
   {
      ::HcalQA *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::HcalQA));
      static ::ROOT::TGenericClassInfo 
         instance("HcalQA", "analyzeRun24or25auau.cpp", 1857,
                  typeid(::HcalQA), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &HcalQA_Dictionary, isa_proxy, 4,
                  sizeof(::HcalQA) );
      instance.SetDelete(&delete_HcalQA);
      instance.SetDeleteArray(&deleteArray_HcalQA);
      instance.SetDestructor(&destruct_HcalQA);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HcalQA*)
   {
      return GenerateInitInstanceLocal(static_cast<::HcalQA*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::HcalQA*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *HcalQA_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::HcalQA*>(nullptr))->GetClass();
      HcalQA_TClassManip(theClass);
   return theClass;
   }

   static void HcalQA_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *SepdPlaneQA_Dictionary();
   static void SepdPlaneQA_TClassManip(TClass*);
   static void delete_SepdPlaneQA(void *p);
   static void deleteArray_SepdPlaneQA(void *p);
   static void destruct_SepdPlaneQA(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SepdPlaneQA*)
   {
      ::SepdPlaneQA *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SepdPlaneQA));
      static ::ROOT::TGenericClassInfo 
         instance("SepdPlaneQA", "analyzeRun24or25auau.cpp", 2033,
                  typeid(::SepdPlaneQA), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &SepdPlaneQA_Dictionary, isa_proxy, 4,
                  sizeof(::SepdPlaneQA) );
      instance.SetDelete(&delete_SepdPlaneQA);
      instance.SetDeleteArray(&deleteArray_SepdPlaneQA);
      instance.SetDestructor(&destruct_SepdPlaneQA);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SepdPlaneQA*)
   {
      return GenerateInitInstanceLocal(static_cast<::SepdPlaneQA*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::SepdPlaneQA*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *SepdPlaneQA_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::SepdPlaneQA*>(nullptr))->GetClass();
      SepdPlaneQA_TClassManip(theClass);
   return theClass;
   }

   static void SepdPlaneQA_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *NSDetectorQAlEMBDTaggR_Dictionary();
   static void NSDetectorQAlEMBDTaggR_TClassManip(TClass*);
   static void delete_NSDetectorQAlEMBDTaggR(void *p);
   static void deleteArray_NSDetectorQAlEMBDTaggR(void *p);
   static void destruct_NSDetectorQAlEMBDTaggR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NSDetectorQA<MBDTag>*)
   {
      ::NSDetectorQA<MBDTag> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::NSDetectorQA<MBDTag>));
      static ::ROOT::TGenericClassInfo 
         instance("NSDetectorQA<MBDTag>", "analyzeRun24or25auau.cpp", 2086,
                  typeid(::NSDetectorQA<MBDTag>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &NSDetectorQAlEMBDTaggR_Dictionary, isa_proxy, 4,
                  sizeof(::NSDetectorQA<MBDTag>) );
      instance.SetDelete(&delete_NSDetectorQAlEMBDTaggR);
      instance.SetDeleteArray(&deleteArray_NSDetectorQAlEMBDTaggR);
      instance.SetDestructor(&destruct_NSDetectorQAlEMBDTaggR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NSDetectorQA<MBDTag>*)
   {
      return GenerateInitInstanceLocal(static_cast<::NSDetectorQA<MBDTag>*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::NSDetectorQA<MBDTag>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *NSDetectorQAlEMBDTaggR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::NSDetectorQA<MBDTag>*>(nullptr))->GetClass();
      NSDetectorQAlEMBDTaggR_TClassManip(theClass);
   return theClass;
   }

   static void NSDetectorQAlEMBDTaggR_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *NSDetectorQAlEsEPDTaggR_Dictionary();
   static void NSDetectorQAlEsEPDTaggR_TClassManip(TClass*);
   static void delete_NSDetectorQAlEsEPDTaggR(void *p);
   static void deleteArray_NSDetectorQAlEsEPDTaggR(void *p);
   static void destruct_NSDetectorQAlEsEPDTaggR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NSDetectorQA<sEPDTag>*)
   {
      ::NSDetectorQA<sEPDTag> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::NSDetectorQA<sEPDTag>));
      static ::ROOT::TGenericClassInfo 
         instance("NSDetectorQA<sEPDTag>", "analyzeRun24or25auau.cpp", 2086,
                  typeid(::NSDetectorQA<sEPDTag>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &NSDetectorQAlEsEPDTaggR_Dictionary, isa_proxy, 4,
                  sizeof(::NSDetectorQA<sEPDTag>) );
      instance.SetDelete(&delete_NSDetectorQAlEsEPDTaggR);
      instance.SetDeleteArray(&deleteArray_NSDetectorQAlEsEPDTaggR);
      instance.SetDestructor(&destruct_NSDetectorQAlEsEPDTaggR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NSDetectorQA<sEPDTag>*)
   {
      return GenerateInitInstanceLocal(static_cast<::NSDetectorQA<sEPDTag>*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::NSDetectorQA<sEPDTag>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *NSDetectorQAlEsEPDTaggR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::NSDetectorQA<sEPDTag>*>(nullptr))->GetClass();
      NSDetectorQAlEsEPDTaggR_TClassManip(theClass);
   return theClass;
   }

   static void NSDetectorQAlEsEPDTaggR_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *MBDTag_Dictionary();
   static void MBDTag_TClassManip(TClass*);
   static void *new_MBDTag(void *p = nullptr);
   static void *newArray_MBDTag(Long_t size, void *p);
   static void delete_MBDTag(void *p);
   static void deleteArray_MBDTag(void *p);
   static void destruct_MBDTag(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MBDTag*)
   {
      ::MBDTag *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::MBDTag));
      static ::ROOT::TGenericClassInfo 
         instance("MBDTag", "analyzeRun24or25auau.cpp", 2206,
                  typeid(::MBDTag), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &MBDTag_Dictionary, isa_proxy, 4,
                  sizeof(::MBDTag) );
      instance.SetNew(&new_MBDTag);
      instance.SetNewArray(&newArray_MBDTag);
      instance.SetDelete(&delete_MBDTag);
      instance.SetDeleteArray(&deleteArray_MBDTag);
      instance.SetDestructor(&destruct_MBDTag);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MBDTag*)
   {
      return GenerateInitInstanceLocal(static_cast<::MBDTag*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::MBDTag*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *MBDTag_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::MBDTag*>(nullptr))->GetClass();
      MBDTag_TClassManip(theClass);
   return theClass;
   }

   static void MBDTag_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *sEPDTag_Dictionary();
   static void sEPDTag_TClassManip(TClass*);
   static void *new_sEPDTag(void *p = nullptr);
   static void *newArray_sEPDTag(Long_t size, void *p);
   static void delete_sEPDTag(void *p);
   static void deleteArray_sEPDTag(void *p);
   static void destruct_sEPDTag(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::sEPDTag*)
   {
      ::sEPDTag *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::sEPDTag));
      static ::ROOT::TGenericClassInfo 
         instance("sEPDTag", "analyzeRun24or25auau.cpp", 2222,
                  typeid(::sEPDTag), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &sEPDTag_Dictionary, isa_proxy, 4,
                  sizeof(::sEPDTag) );
      instance.SetNew(&new_sEPDTag);
      instance.SetNewArray(&newArray_sEPDTag);
      instance.SetDelete(&delete_sEPDTag);
      instance.SetDeleteArray(&deleteArray_sEPDTag);
      instance.SetDestructor(&destruct_sEPDTag);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::sEPDTag*)
   {
      return GenerateInitInstanceLocal(static_cast<::sEPDTag*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::sEPDTag*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *sEPDTag_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::sEPDTag*>(nullptr))->GetClass();
      sEPDTag_TClassManip(theClass);
   return theClass;
   }

   static void sEPDTag_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *EventQA_Dictionary();
   static void EventQA_TClassManip(TClass*);
   static void delete_EventQA(void *p);
   static void deleteArray_EventQA(void *p);
   static void destruct_EventQA(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EventQA*)
   {
      ::EventQA *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::EventQA));
      static ::ROOT::TGenericClassInfo 
         instance("EventQA", "analyzeRun24or25auau.cpp", 2247,
                  typeid(::EventQA), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &EventQA_Dictionary, isa_proxy, 4,
                  sizeof(::EventQA) );
      instance.SetDelete(&delete_EventQA);
      instance.SetDeleteArray(&deleteArray_EventQA);
      instance.SetDestructor(&destruct_EventQA);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EventQA*)
   {
      return GenerateInitInstanceLocal(static_cast<::EventQA*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::EventQA*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *EventQA_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::EventQA*>(nullptr))->GetClass();
      EventQA_TClassManip(theClass);
   return theClass;
   }

   static void EventQA_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *EventQAcLcLVzPoint_Dictionary();
   static void EventQAcLcLVzPoint_TClassManip(TClass*);
   static void *new_EventQAcLcLVzPoint(void *p = nullptr);
   static void *newArray_EventQAcLcLVzPoint(Long_t size, void *p);
   static void delete_EventQAcLcLVzPoint(void *p);
   static void deleteArray_EventQAcLcLVzPoint(void *p);
   static void destruct_EventQAcLcLVzPoint(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EventQA::VzPoint*)
   {
      ::EventQA::VzPoint *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::EventQA::VzPoint));
      static ::ROOT::TGenericClassInfo 
         instance("EventQA::VzPoint", "analyzeRun24or25auau.cpp", 2503,
                  typeid(::EventQA::VzPoint), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &EventQAcLcLVzPoint_Dictionary, isa_proxy, 4,
                  sizeof(::EventQA::VzPoint) );
      instance.SetNew(&new_EventQAcLcLVzPoint);
      instance.SetNewArray(&newArray_EventQAcLcLVzPoint);
      instance.SetDelete(&delete_EventQAcLcLVzPoint);
      instance.SetDeleteArray(&deleteArray_EventQAcLcLVzPoint);
      instance.SetDestructor(&destruct_EventQAcLcLVzPoint);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EventQA::VzPoint*)
   {
      return GenerateInitInstanceLocal(static_cast<::EventQA::VzPoint*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::EventQA::VzPoint*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *EventQAcLcLVzPoint_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::EventQA::VzPoint*>(nullptr))->GetClass();
      EventQAcLcLVzPoint_TClassManip(theClass);
   return theClass;
   }

   static void EventQAcLcLVzPoint_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *JetQA_Dictionary();
   static void JetQA_TClassManip(TClass*);
   static void delete_JetQA(void *p);
   static void deleteArray_JetQA(void *p);
   static void destruct_JetQA(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::JetQA*)
   {
      ::JetQA *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::JetQA));
      static ::ROOT::TGenericClassInfo 
         instance("JetQA", "analyzeRun24or25auau.cpp", 2528,
                  typeid(::JetQA), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &JetQA_Dictionary, isa_proxy, 4,
                  sizeof(::JetQA) );
      instance.SetDelete(&delete_JetQA);
      instance.SetDeleteArray(&deleteArray_JetQA);
      instance.SetDestructor(&destruct_JetQA);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::JetQA*)
   {
      return GenerateInitInstanceLocal(static_cast<::JetQA*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::JetQA*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *JetQA_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::JetQA*>(nullptr))->GetClass();
      JetQA_TClassManip(theClass);
   return theClass;
   }

   static void JetQA_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *VnPlotQA_Dictionary();
   static void VnPlotQA_TClassManip(TClass*);
   static void delete_VnPlotQA(void *p);
   static void deleteArray_VnPlotQA(void *p);
   static void destruct_VnPlotQA(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::VnPlotQA*)
   {
      ::VnPlotQA *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::VnPlotQA));
      static ::ROOT::TGenericClassInfo 
         instance("VnPlotQA", "analyzeRun24or25auau.cpp", 2767,
                  typeid(::VnPlotQA), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &VnPlotQA_Dictionary, isa_proxy, 4,
                  sizeof(::VnPlotQA) );
      instance.SetDelete(&delete_VnPlotQA);
      instance.SetDeleteArray(&deleteArray_VnPlotQA);
      instance.SetDestructor(&destruct_VnPlotQA);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::VnPlotQA*)
   {
      return GenerateInitInstanceLocal(static_cast<::VnPlotQA*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::VnPlotQA*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *VnPlotQA_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::VnPlotQA*>(nullptr))->GetClass();
      VnPlotQA_TClassManip(theClass);
   return theClass;
   }

   static void VnPlotQA_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/./analyzeRun24or25auau.cpp");
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_CutKey(void *p) {
      return  p ? new(p) ::CutKey : new ::CutKey;
   }
   static void *newArray_CutKey(Long_t nElements, void *p) {
      return p ? new(p) ::CutKey[nElements] : new ::CutKey[nElements];
   }
   // Wrapper around operator delete
   static void delete_CutKey(void *p) {
      delete (static_cast<::CutKey*>(p));
   }
   static void deleteArray_CutKey(void *p) {
      delete [] (static_cast<::CutKey*>(p));
   }
   static void destruct_CutKey(void *p) {
      typedef ::CutKey current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::CutKey

namespace ROOT {
   // Wrappers around operator new
   static void *new_NSCachelEMapPairgR(void *p) {
      return  p ? new(p) ::NSCache<MapPair> : new ::NSCache<MapPair>;
   }
   static void *newArray_NSCachelEMapPairgR(Long_t nElements, void *p) {
      return p ? new(p) ::NSCache<MapPair>[nElements] : new ::NSCache<MapPair>[nElements];
   }
   // Wrapper around operator delete
   static void delete_NSCachelEMapPairgR(void *p) {
      delete (static_cast<::NSCache<MapPair>*>(p));
   }
   static void deleteArray_NSCachelEMapPairgR(void *p) {
      delete [] (static_cast<::NSCache<MapPair>*>(p));
   }
   static void destruct_NSCachelEMapPairgR(void *p) {
      typedef ::NSCache<MapPair> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::NSCache<MapPair>

namespace ROOT {
   // Wrappers around operator new
   static void *new_NSCachelEMapPairgRcLcLH(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::NSCache<MapPair>::H : new ::NSCache<MapPair>::H;
   }
   static void *newArray_NSCachelEMapPairgRcLcLH(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::NSCache<MapPair>::H[nElements] : new ::NSCache<MapPair>::H[nElements];
   }
   // Wrapper around operator delete
   static void delete_NSCachelEMapPairgRcLcLH(void *p) {
      delete (static_cast<::NSCache<MapPair>::H*>(p));
   }
   static void deleteArray_NSCachelEMapPairgRcLcLH(void *p) {
      delete [] (static_cast<::NSCache<MapPair>::H*>(p));
   }
   static void destruct_NSCachelEMapPairgRcLcLH(void *p) {
      typedef ::NSCache<MapPair>::H current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::NSCache<MapPair>::H

namespace ROOT {
   // Wrappers around operator new
   static void *new_MapPair(void *p) {
      return  p ? new(p) ::MapPair : new ::MapPair;
   }
   static void *newArray_MapPair(Long_t nElements, void *p) {
      return p ? new(p) ::MapPair[nElements] : new ::MapPair[nElements];
   }
   // Wrapper around operator delete
   static void delete_MapPair(void *p) {
      delete (static_cast<::MapPair*>(p));
   }
   static void deleteArray_MapPair(void *p) {
      delete [] (static_cast<::MapPair*>(p));
   }
   static void destruct_MapPair(void *p) {
      typedef ::MapPair current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::MapPair

namespace ROOT {
   // Wrapper around operator delete
   static void delete_QA(void *p) {
      delete (static_cast<::QA*>(p));
   }
   static void deleteArray_QA(void *p) {
      delete [] (static_cast<::QA*>(p));
   }
   static void destruct_QA(void *p) {
      typedef ::QA current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::QA

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pi0QA(void *p) {
      delete (static_cast<::Pi0QA*>(p));
   }
   static void deleteArray_Pi0QA(void *p) {
      delete [] (static_cast<::Pi0QA*>(p));
   }
   static void destruct_Pi0QA(void *p) {
      typedef ::Pi0QA current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pi0QA

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pi0QAcLcLFitInfo(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pi0QA::FitInfo : new ::Pi0QA::FitInfo;
   }
   static void *newArray_Pi0QAcLcLFitInfo(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pi0QA::FitInfo[nElements] : new ::Pi0QA::FitInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pi0QAcLcLFitInfo(void *p) {
      delete (static_cast<::Pi0QA::FitInfo*>(p));
   }
   static void deleteArray_Pi0QAcLcLFitInfo(void *p) {
      delete [] (static_cast<::Pi0QA::FitInfo*>(p));
   }
   static void destruct_Pi0QAcLcLFitInfo(void *p) {
      typedef ::Pi0QA::FitInfo current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pi0QA::FitInfo

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pi0QAcLcLFitPair(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pi0QA::FitPair : new ::Pi0QA::FitPair;
   }
   static void *newArray_Pi0QAcLcLFitPair(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pi0QA::FitPair[nElements] : new ::Pi0QA::FitPair[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pi0QAcLcLFitPair(void *p) {
      delete (static_cast<::Pi0QA::FitPair*>(p));
   }
   static void deleteArray_Pi0QAcLcLFitPair(void *p) {
      delete [] (static_cast<::Pi0QA::FitPair*>(p));
   }
   static void destruct_Pi0QAcLcLFitPair(void *p) {
      typedef ::Pi0QA::FitPair current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pi0QA::FitPair

namespace ROOT {
   // Wrappers around operator new
   static void *new_NSPair(void *p) {
      return  p ? new(p) ::NSPair : new ::NSPair;
   }
   static void *newArray_NSPair(Long_t nElements, void *p) {
      return p ? new(p) ::NSPair[nElements] : new ::NSPair[nElements];
   }
   // Wrapper around operator delete
   static void delete_NSPair(void *p) {
      delete (static_cast<::NSPair*>(p));
   }
   static void deleteArray_NSPair(void *p) {
      delete [] (static_cast<::NSPair*>(p));
   }
   static void destruct_NSPair(void *p) {
      typedef ::NSPair current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::NSPair

namespace ROOT {
   // Wrapper around operator delete
   static void delete_CorrQA(void *p) {
      delete (static_cast<::CorrQA*>(p));
   }
   static void deleteArray_CorrQA(void *p) {
      delete [] (static_cast<::CorrQA*>(p));
   }
   static void destruct_CorrQA(void *p) {
      typedef ::CorrQA current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::CorrQA

namespace ROOT {
   // Wrapper around operator delete
   static void delete_EmcalQA(void *p) {
      delete (static_cast<::EmcalQA*>(p));
   }
   static void deleteArray_EmcalQA(void *p) {
      delete [] (static_cast<::EmcalQA*>(p));
   }
   static void destruct_EmcalQA(void *p) {
      typedef ::EmcalQA current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::EmcalQA

namespace ROOT {
   // Wrapper around operator delete
   static void delete_HcalQA(void *p) {
      delete (static_cast<::HcalQA*>(p));
   }
   static void deleteArray_HcalQA(void *p) {
      delete [] (static_cast<::HcalQA*>(p));
   }
   static void destruct_HcalQA(void *p) {
      typedef ::HcalQA current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::HcalQA

namespace ROOT {
   // Wrapper around operator delete
   static void delete_SepdPlaneQA(void *p) {
      delete (static_cast<::SepdPlaneQA*>(p));
   }
   static void deleteArray_SepdPlaneQA(void *p) {
      delete [] (static_cast<::SepdPlaneQA*>(p));
   }
   static void destruct_SepdPlaneQA(void *p) {
      typedef ::SepdPlaneQA current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::SepdPlaneQA

namespace ROOT {
   // Wrapper around operator delete
   static void delete_NSDetectorQAlEMBDTaggR(void *p) {
      delete (static_cast<::NSDetectorQA<MBDTag>*>(p));
   }
   static void deleteArray_NSDetectorQAlEMBDTaggR(void *p) {
      delete [] (static_cast<::NSDetectorQA<MBDTag>*>(p));
   }
   static void destruct_NSDetectorQAlEMBDTaggR(void *p) {
      typedef ::NSDetectorQA<MBDTag> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::NSDetectorQA<MBDTag>

namespace ROOT {
   // Wrapper around operator delete
   static void delete_NSDetectorQAlEsEPDTaggR(void *p) {
      delete (static_cast<::NSDetectorQA<sEPDTag>*>(p));
   }
   static void deleteArray_NSDetectorQAlEsEPDTaggR(void *p) {
      delete [] (static_cast<::NSDetectorQA<sEPDTag>*>(p));
   }
   static void destruct_NSDetectorQAlEsEPDTaggR(void *p) {
      typedef ::NSDetectorQA<sEPDTag> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::NSDetectorQA<sEPDTag>

namespace ROOT {
   // Wrappers around operator new
   static void *new_MBDTag(void *p) {
      return  p ? new(p) ::MBDTag : new ::MBDTag;
   }
   static void *newArray_MBDTag(Long_t nElements, void *p) {
      return p ? new(p) ::MBDTag[nElements] : new ::MBDTag[nElements];
   }
   // Wrapper around operator delete
   static void delete_MBDTag(void *p) {
      delete (static_cast<::MBDTag*>(p));
   }
   static void deleteArray_MBDTag(void *p) {
      delete [] (static_cast<::MBDTag*>(p));
   }
   static void destruct_MBDTag(void *p) {
      typedef ::MBDTag current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::MBDTag

namespace ROOT {
   // Wrappers around operator new
   static void *new_sEPDTag(void *p) {
      return  p ? new(p) ::sEPDTag : new ::sEPDTag;
   }
   static void *newArray_sEPDTag(Long_t nElements, void *p) {
      return p ? new(p) ::sEPDTag[nElements] : new ::sEPDTag[nElements];
   }
   // Wrapper around operator delete
   static void delete_sEPDTag(void *p) {
      delete (static_cast<::sEPDTag*>(p));
   }
   static void deleteArray_sEPDTag(void *p) {
      delete [] (static_cast<::sEPDTag*>(p));
   }
   static void destruct_sEPDTag(void *p) {
      typedef ::sEPDTag current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::sEPDTag

namespace ROOT {
   // Wrapper around operator delete
   static void delete_EventQA(void *p) {
      delete (static_cast<::EventQA*>(p));
   }
   static void deleteArray_EventQA(void *p) {
      delete [] (static_cast<::EventQA*>(p));
   }
   static void destruct_EventQA(void *p) {
      typedef ::EventQA current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::EventQA

namespace ROOT {
   // Wrappers around operator new
   static void *new_EventQAcLcLVzPoint(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::EventQA::VzPoint : new ::EventQA::VzPoint;
   }
   static void *newArray_EventQAcLcLVzPoint(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::EventQA::VzPoint[nElements] : new ::EventQA::VzPoint[nElements];
   }
   // Wrapper around operator delete
   static void delete_EventQAcLcLVzPoint(void *p) {
      delete (static_cast<::EventQA::VzPoint*>(p));
   }
   static void deleteArray_EventQAcLcLVzPoint(void *p) {
      delete [] (static_cast<::EventQA::VzPoint*>(p));
   }
   static void destruct_EventQAcLcLVzPoint(void *p) {
      typedef ::EventQA::VzPoint current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::EventQA::VzPoint

namespace ROOT {
   // Wrapper around operator delete
   static void delete_JetQA(void *p) {
      delete (static_cast<::JetQA*>(p));
   }
   static void deleteArray_JetQA(void *p) {
      delete [] (static_cast<::JetQA*>(p));
   }
   static void destruct_JetQA(void *p) {
      typedef ::JetQA current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::JetQA

namespace ROOT {
   // Wrapper around operator delete
   static void delete_VnPlotQA(void *p) {
      delete (static_cast<::VnPlotQA*>(p));
   }
   static void deleteArray_VnPlotQA(void *p) {
      delete [] (static_cast<::VnPlotQA*>(p));
   }
   static void destruct_VnPlotQA(void *p) {
      typedef ::VnPlotQA current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::VnPlotQA

namespace ROOT {
   static TClass *vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR_Dictionary();
   static void vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR(Long_t size, void *p);
   static void delete_vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR(void *p);
   static void deleteArray_vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR(void *p);
   static void destruct_vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<pair<string,shared_ptr<TH2> > >*)
   {
      vector<pair<string,shared_ptr<TH2> > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<pair<string,shared_ptr<TH2> > >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<pair<string,shared_ptr<TH2> > >", -2, "vector", 387,
                  typeid(vector<pair<string,shared_ptr<TH2> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<pair<string,shared_ptr<TH2> > >) );
      instance.SetNew(&new_vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR);
      instance.SetNewArray(&newArray_vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR);
      instance.SetDelete(&delete_vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR);
      instance.SetDestructor(&destruct_vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<pair<string,shared_ptr<TH2> > > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<pair<string,shared_ptr<TH2> > >","std::__1::vector<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::shared_ptr<TH2>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::shared_ptr<TH2>>>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<pair<string,shared_ptr<TH2> > >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<pair<string,shared_ptr<TH2> > >*>(nullptr))->GetClass();
      vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<pair<string,shared_ptr<TH2> > > : new vector<pair<string,shared_ptr<TH2> > >;
   }
   static void *newArray_vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<pair<string,shared_ptr<TH2> > >[nElements] : new vector<pair<string,shared_ptr<TH2> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR(void *p) {
      delete (static_cast<vector<pair<string,shared_ptr<TH2> > >*>(p));
   }
   static void deleteArray_vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR(void *p) {
      delete [] (static_cast<vector<pair<string,shared_ptr<TH2> > >*>(p));
   }
   static void destruct_vectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgR(void *p) {
      typedef vector<pair<string,shared_ptr<TH2> > > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<pair<string,shared_ptr<TH2> > >

namespace ROOT {
   static TClass *vectorlEpairlEfloatcOfloatgRsPgR_Dictionary();
   static void vectorlEpairlEfloatcOfloatgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEpairlEfloatcOfloatgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEpairlEfloatcOfloatgRsPgR(Long_t size, void *p);
   static void delete_vectorlEpairlEfloatcOfloatgRsPgR(void *p);
   static void deleteArray_vectorlEpairlEfloatcOfloatgRsPgR(void *p);
   static void destruct_vectorlEpairlEfloatcOfloatgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<pair<float,float> >*)
   {
      vector<pair<float,float> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<pair<float,float> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<pair<float,float> >", -2, "vector", 387,
                  typeid(vector<pair<float,float> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEpairlEfloatcOfloatgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<pair<float,float> >) );
      instance.SetNew(&new_vectorlEpairlEfloatcOfloatgRsPgR);
      instance.SetNewArray(&newArray_vectorlEpairlEfloatcOfloatgRsPgR);
      instance.SetDelete(&delete_vectorlEpairlEfloatcOfloatgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEpairlEfloatcOfloatgRsPgR);
      instance.SetDestructor(&destruct_vectorlEpairlEfloatcOfloatgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<pair<float,float> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<pair<float,float> >","std::__1::vector<std::__1::pair<float, float>, std::__1::allocator<std::__1::pair<float, float>>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<pair<float,float> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEpairlEfloatcOfloatgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<pair<float,float> >*>(nullptr))->GetClass();
      vectorlEpairlEfloatcOfloatgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEpairlEfloatcOfloatgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEpairlEfloatcOfloatgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<pair<float,float> > : new vector<pair<float,float> >;
   }
   static void *newArray_vectorlEpairlEfloatcOfloatgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<pair<float,float> >[nElements] : new vector<pair<float,float> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEpairlEfloatcOfloatgRsPgR(void *p) {
      delete (static_cast<vector<pair<float,float> >*>(p));
   }
   static void deleteArray_vectorlEpairlEfloatcOfloatgRsPgR(void *p) {
      delete [] (static_cast<vector<pair<float,float> >*>(p));
   }
   static void destruct_vectorlEpairlEfloatcOfloatgRsPgR(void *p) {
      typedef vector<pair<float,float> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<pair<float,float> >

namespace ROOT {
   static TClass *vectorlETProfilemUgR_Dictionary();
   static void vectorlETProfilemUgR_TClassManip(TClass*);
   static void *new_vectorlETProfilemUgR(void *p = nullptr);
   static void *newArray_vectorlETProfilemUgR(Long_t size, void *p);
   static void delete_vectorlETProfilemUgR(void *p);
   static void deleteArray_vectorlETProfilemUgR(void *p);
   static void destruct_vectorlETProfilemUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TProfile*>*)
   {
      vector<TProfile*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TProfile*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TProfile*>", -2, "vector", 387,
                  typeid(vector<TProfile*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETProfilemUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TProfile*>) );
      instance.SetNew(&new_vectorlETProfilemUgR);
      instance.SetNewArray(&newArray_vectorlETProfilemUgR);
      instance.SetDelete(&delete_vectorlETProfilemUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETProfilemUgR);
      instance.SetDestructor(&destruct_vectorlETProfilemUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TProfile*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<TProfile*>","std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TProfile*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETProfilemUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TProfile*>*>(nullptr))->GetClass();
      vectorlETProfilemUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETProfilemUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETProfilemUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TProfile*> : new vector<TProfile*>;
   }
   static void *newArray_vectorlETProfilemUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TProfile*>[nElements] : new vector<TProfile*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETProfilemUgR(void *p) {
      delete (static_cast<vector<TProfile*>*>(p));
   }
   static void deleteArray_vectorlETProfilemUgR(void *p) {
      delete [] (static_cast<vector<TProfile*>*>(p));
   }
   static void destruct_vectorlETProfilemUgR(void *p) {
      typedef vector<TProfile*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TProfile*>

namespace ROOT {
   static TClass *vectorlETH1mUgR_Dictionary();
   static void vectorlETH1mUgR_TClassManip(TClass*);
   static void *new_vectorlETH1mUgR(void *p = nullptr);
   static void *newArray_vectorlETH1mUgR(Long_t size, void *p);
   static void delete_vectorlETH1mUgR(void *p);
   static void deleteArray_vectorlETH1mUgR(void *p);
   static void destruct_vectorlETH1mUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TH1*>*)
   {
      vector<TH1*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TH1*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TH1*>", -2, "vector", 387,
                  typeid(vector<TH1*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETH1mUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TH1*>) );
      instance.SetNew(&new_vectorlETH1mUgR);
      instance.SetNewArray(&newArray_vectorlETH1mUgR);
      instance.SetDelete(&delete_vectorlETH1mUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETH1mUgR);
      instance.SetDestructor(&destruct_vectorlETH1mUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TH1*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<TH1*>","std::__1::vector<TH1*, std::__1::allocator<TH1*>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TH1*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETH1mUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TH1*>*>(nullptr))->GetClass();
      vectorlETH1mUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETH1mUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETH1mUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TH1*> : new vector<TH1*>;
   }
   static void *newArray_vectorlETH1mUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TH1*>[nElements] : new vector<TH1*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETH1mUgR(void *p) {
      delete (static_cast<vector<TH1*>*>(p));
   }
   static void deleteArray_vectorlETH1mUgR(void *p) {
      delete [] (static_cast<vector<TH1*>*>(p));
   }
   static void destruct_vectorlETH1mUgR(void *p) {
      typedef vector<TH1*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TH1*>

namespace ROOT {
   static TClass *vectorlETGraphErrorsmUgR_Dictionary();
   static void vectorlETGraphErrorsmUgR_TClassManip(TClass*);
   static void *new_vectorlETGraphErrorsmUgR(void *p = nullptr);
   static void *newArray_vectorlETGraphErrorsmUgR(Long_t size, void *p);
   static void delete_vectorlETGraphErrorsmUgR(void *p);
   static void deleteArray_vectorlETGraphErrorsmUgR(void *p);
   static void destruct_vectorlETGraphErrorsmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TGraphErrors*>*)
   {
      vector<TGraphErrors*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TGraphErrors*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TGraphErrors*>", -2, "vector", 387,
                  typeid(vector<TGraphErrors*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETGraphErrorsmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TGraphErrors*>) );
      instance.SetNew(&new_vectorlETGraphErrorsmUgR);
      instance.SetNewArray(&newArray_vectorlETGraphErrorsmUgR);
      instance.SetDelete(&delete_vectorlETGraphErrorsmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETGraphErrorsmUgR);
      instance.SetDestructor(&destruct_vectorlETGraphErrorsmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TGraphErrors*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<TGraphErrors*>","std::__1::vector<TGraphErrors*, std::__1::allocator<TGraphErrors*>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TGraphErrors*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETGraphErrorsmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TGraphErrors*>*>(nullptr))->GetClass();
      vectorlETGraphErrorsmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETGraphErrorsmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETGraphErrorsmUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TGraphErrors*> : new vector<TGraphErrors*>;
   }
   static void *newArray_vectorlETGraphErrorsmUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TGraphErrors*>[nElements] : new vector<TGraphErrors*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETGraphErrorsmUgR(void *p) {
      delete (static_cast<vector<TGraphErrors*>*>(p));
   }
   static void deleteArray_vectorlETGraphErrorsmUgR(void *p) {
      delete [] (static_cast<vector<TGraphErrors*>*>(p));
   }
   static void destruct_vectorlETGraphErrorsmUgR(void *p) {
      typedef vector<TGraphErrors*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TGraphErrors*>

namespace ROOT {
   static TClass *unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR_Dictionary();
   static void unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR_TClassManip(TClass*);
   static void *new_unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR(void *p = nullptr);
   static void *newArray_unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR(Long_t size, void *p);
   static void delete_unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR(void *p);
   static void deleteArray_unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR(void *p);
   static void destruct_unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const unordered_map<string,vector<pair<string,shared_ptr<TH2> > > >*)
   {
      unordered_map<string,vector<pair<string,shared_ptr<TH2> > > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(unordered_map<string,vector<pair<string,shared_ptr<TH2> > > >));
      static ::ROOT::TGenericClassInfo 
         instance("unordered_map<string,vector<pair<string,shared_ptr<TH2> > > >", -2, "unordered_map", 1029,
                  typeid(unordered_map<string,vector<pair<string,shared_ptr<TH2> > > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(unordered_map<string,vector<pair<string,shared_ptr<TH2> > > >) );
      instance.SetNew(&new_unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR);
      instance.SetNewArray(&newArray_unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR);
      instance.SetDelete(&delete_unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR);
      instance.SetDestructor(&destruct_unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< unordered_map<string,vector<pair<string,shared_ptr<TH2> > > > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("unordered_map<string,vector<pair<string,shared_ptr<TH2> > > >","std::__1::unordered_map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::vector<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::shared_ptr<TH2>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::shared_ptr<TH2>>>>, std::__1::hash<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::equal_to<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::vector<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::shared_ptr<TH2>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::shared_ptr<TH2>>>>>>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const unordered_map<string,vector<pair<string,shared_ptr<TH2> > > >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const unordered_map<string,vector<pair<string,shared_ptr<TH2> > > >*>(nullptr))->GetClass();
      unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,vector<pair<string,shared_ptr<TH2> > > > : new unordered_map<string,vector<pair<string,shared_ptr<TH2> > > >;
   }
   static void *newArray_unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,vector<pair<string,shared_ptr<TH2> > > >[nElements] : new unordered_map<string,vector<pair<string,shared_ptr<TH2> > > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR(void *p) {
      delete (static_cast<unordered_map<string,vector<pair<string,shared_ptr<TH2> > > >*>(p));
   }
   static void deleteArray_unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR(void *p) {
      delete [] (static_cast<unordered_map<string,vector<pair<string,shared_ptr<TH2> > > >*>(p));
   }
   static void destruct_unordered_maplEstringcOvectorlEpairlEstringcOshared_ptrlETH2gRsPgRsPgRsPgR(void *p) {
      typedef unordered_map<string,vector<pair<string,shared_ptr<TH2> > > > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class unordered_map<string,vector<pair<string,shared_ptr<TH2> > > >

namespace ROOT {
   static TClass *unordered_maplEstringcOvectorlETH1mUgRsPgR_Dictionary();
   static void unordered_maplEstringcOvectorlETH1mUgRsPgR_TClassManip(TClass*);
   static void *new_unordered_maplEstringcOvectorlETH1mUgRsPgR(void *p = nullptr);
   static void *newArray_unordered_maplEstringcOvectorlETH1mUgRsPgR(Long_t size, void *p);
   static void delete_unordered_maplEstringcOvectorlETH1mUgRsPgR(void *p);
   static void deleteArray_unordered_maplEstringcOvectorlETH1mUgRsPgR(void *p);
   static void destruct_unordered_maplEstringcOvectorlETH1mUgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const unordered_map<string,vector<TH1*> >*)
   {
      unordered_map<string,vector<TH1*> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(unordered_map<string,vector<TH1*> >));
      static ::ROOT::TGenericClassInfo 
         instance("unordered_map<string,vector<TH1*> >", -2, "unordered_map", 1029,
                  typeid(unordered_map<string,vector<TH1*> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &unordered_maplEstringcOvectorlETH1mUgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(unordered_map<string,vector<TH1*> >) );
      instance.SetNew(&new_unordered_maplEstringcOvectorlETH1mUgRsPgR);
      instance.SetNewArray(&newArray_unordered_maplEstringcOvectorlETH1mUgRsPgR);
      instance.SetDelete(&delete_unordered_maplEstringcOvectorlETH1mUgRsPgR);
      instance.SetDeleteArray(&deleteArray_unordered_maplEstringcOvectorlETH1mUgRsPgR);
      instance.SetDestructor(&destruct_unordered_maplEstringcOvectorlETH1mUgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< unordered_map<string,vector<TH1*> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("unordered_map<string,vector<TH1*> >","std::__1::unordered_map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::vector<TH1*, std::__1::allocator<TH1*>>, std::__1::hash<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::equal_to<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::vector<TH1*, std::__1::allocator<TH1*>>>>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const unordered_map<string,vector<TH1*> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *unordered_maplEstringcOvectorlETH1mUgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const unordered_map<string,vector<TH1*> >*>(nullptr))->GetClass();
      unordered_maplEstringcOvectorlETH1mUgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void unordered_maplEstringcOvectorlETH1mUgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_unordered_maplEstringcOvectorlETH1mUgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,vector<TH1*> > : new unordered_map<string,vector<TH1*> >;
   }
   static void *newArray_unordered_maplEstringcOvectorlETH1mUgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,vector<TH1*> >[nElements] : new unordered_map<string,vector<TH1*> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_unordered_maplEstringcOvectorlETH1mUgRsPgR(void *p) {
      delete (static_cast<unordered_map<string,vector<TH1*> >*>(p));
   }
   static void deleteArray_unordered_maplEstringcOvectorlETH1mUgRsPgR(void *p) {
      delete [] (static_cast<unordered_map<string,vector<TH1*> >*>(p));
   }
   static void destruct_unordered_maplEstringcOvectorlETH1mUgRsPgR(void *p) {
      typedef unordered_map<string,vector<TH1*> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class unordered_map<string,vector<TH1*> >

namespace ROOT {
   static TClass *unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR_Dictionary();
   static void unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR_TClassManip(TClass*);
   static void *new_unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR(void *p = nullptr);
   static void *newArray_unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR(Long_t size, void *p);
   static void delete_unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR(void *p);
   static void deleteArray_unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR(void *p);
   static void destruct_unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const unordered_map<string,unordered_map<string,map<int,map<string,vector<TProfile*> > > > >*)
   {
      unordered_map<string,unordered_map<string,map<int,map<string,vector<TProfile*> > > > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(unordered_map<string,unordered_map<string,map<int,map<string,vector<TProfile*> > > > >));
      static ::ROOT::TGenericClassInfo 
         instance("unordered_map<string,unordered_map<string,map<int,map<string,vector<TProfile*> > > > >", -2, "unordered_map", 1029,
                  typeid(unordered_map<string,unordered_map<string,map<int,map<string,vector<TProfile*> > > > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(unordered_map<string,unordered_map<string,map<int,map<string,vector<TProfile*> > > > >) );
      instance.SetNew(&new_unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR);
      instance.SetNewArray(&newArray_unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR);
      instance.SetDelete(&delete_unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR);
      instance.SetDestructor(&destruct_unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< unordered_map<string,unordered_map<string,map<int,map<string,vector<TProfile*> > > > > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("unordered_map<string,unordered_map<string,map<int,map<string,vector<TProfile*> > > > >","std::__1::unordered_map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::unordered_map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::map<int, std::__1::map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>, std::__1::less<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>>>>, std::__1::less<int>, std::__1::allocator<std::__1::pair<int const, std::__1::map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>, std::__1::less<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>>>>>>>, std::__1::hash<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::equal_to<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::map<int, std::__1::map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>, std::__1::less<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>>>>, std::__1::less<int>, std::__1::allocator<std::__1::pair<int const, std::__1::map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>, std::__1::less<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>>>>>>>>>>, std::__1::hash<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::equal_to<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::unordered_map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::map<int, std::__1::map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>, std::__1::less<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>>>>, std::__1::less<int>, std::__1::allocator<std::__1::pair<int const, std::__1::map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>, std::__1::less<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>>>>>>>, std::__1::hash<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::equal_to<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::map<int, std::__1::map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>, std::__1::less<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>>>>, std::__1::less<int>, std::__1::allocator<std::__1::pair<int const, std::__1::map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>, std::__1::less<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>>>>>>>>>>>>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const unordered_map<string,unordered_map<string,map<int,map<string,vector<TProfile*> > > > >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const unordered_map<string,unordered_map<string,map<int,map<string,vector<TProfile*> > > > >*>(nullptr))->GetClass();
      unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,unordered_map<string,map<int,map<string,vector<TProfile*> > > > > : new unordered_map<string,unordered_map<string,map<int,map<string,vector<TProfile*> > > > >;
   }
   static void *newArray_unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,unordered_map<string,map<int,map<string,vector<TProfile*> > > > >[nElements] : new unordered_map<string,unordered_map<string,map<int,map<string,vector<TProfile*> > > > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR(void *p) {
      delete (static_cast<unordered_map<string,unordered_map<string,map<int,map<string,vector<TProfile*> > > > >*>(p));
   }
   static void deleteArray_unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR(void *p) {
      delete [] (static_cast<unordered_map<string,unordered_map<string,map<int,map<string,vector<TProfile*> > > > >*>(p));
   }
   static void destruct_unordered_maplEstringcOunordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgRsPgR(void *p) {
      typedef unordered_map<string,unordered_map<string,map<int,map<string,vector<TProfile*> > > > > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class unordered_map<string,unordered_map<string,map<int,map<string,vector<TProfile*> > > > >

namespace ROOT {
   static TClass *unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR_Dictionary();
   static void unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR_TClassManip(TClass*);
   static void *new_unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR(void *p = nullptr);
   static void *newArray_unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR(Long_t size, void *p);
   static void delete_unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR(void *p);
   static void deleteArray_unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR(void *p);
   static void destruct_unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const unordered_map<string,map<int,map<string,vector<TProfile*> > > >*)
   {
      unordered_map<string,map<int,map<string,vector<TProfile*> > > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(unordered_map<string,map<int,map<string,vector<TProfile*> > > >));
      static ::ROOT::TGenericClassInfo 
         instance("unordered_map<string,map<int,map<string,vector<TProfile*> > > >", -2, "unordered_map", 1029,
                  typeid(unordered_map<string,map<int,map<string,vector<TProfile*> > > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(unordered_map<string,map<int,map<string,vector<TProfile*> > > >) );
      instance.SetNew(&new_unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR);
      instance.SetNewArray(&newArray_unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR);
      instance.SetDelete(&delete_unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR);
      instance.SetDestructor(&destruct_unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< unordered_map<string,map<int,map<string,vector<TProfile*> > > > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("unordered_map<string,map<int,map<string,vector<TProfile*> > > >","std::__1::unordered_map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::map<int, std::__1::map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>, std::__1::less<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>>>>, std::__1::less<int>, std::__1::allocator<std::__1::pair<int const, std::__1::map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>, std::__1::less<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>>>>>>>, std::__1::hash<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::equal_to<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::map<int, std::__1::map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>, std::__1::less<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>>>>, std::__1::less<int>, std::__1::allocator<std::__1::pair<int const, std::__1::map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>, std::__1::less<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>>>>>>>>>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const unordered_map<string,map<int,map<string,vector<TProfile*> > > >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const unordered_map<string,map<int,map<string,vector<TProfile*> > > >*>(nullptr))->GetClass();
      unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,map<int,map<string,vector<TProfile*> > > > : new unordered_map<string,map<int,map<string,vector<TProfile*> > > >;
   }
   static void *newArray_unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,map<int,map<string,vector<TProfile*> > > >[nElements] : new unordered_map<string,map<int,map<string,vector<TProfile*> > > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR(void *p) {
      delete (static_cast<unordered_map<string,map<int,map<string,vector<TProfile*> > > >*>(p));
   }
   static void deleteArray_unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR(void *p) {
      delete [] (static_cast<unordered_map<string,map<int,map<string,vector<TProfile*> > > >*>(p));
   }
   static void destruct_unordered_maplEstringcOmaplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgRsPgR(void *p) {
      typedef unordered_map<string,map<int,map<string,vector<TProfile*> > > > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class unordered_map<string,map<int,map<string,vector<TProfile*> > > >

namespace ROOT {
   static TClass *unordered_maplEstringcOTH2FmUgR_Dictionary();
   static void unordered_maplEstringcOTH2FmUgR_TClassManip(TClass*);
   static void *new_unordered_maplEstringcOTH2FmUgR(void *p = nullptr);
   static void *newArray_unordered_maplEstringcOTH2FmUgR(Long_t size, void *p);
   static void delete_unordered_maplEstringcOTH2FmUgR(void *p);
   static void deleteArray_unordered_maplEstringcOTH2FmUgR(void *p);
   static void destruct_unordered_maplEstringcOTH2FmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const unordered_map<string,TH2F*>*)
   {
      unordered_map<string,TH2F*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(unordered_map<string,TH2F*>));
      static ::ROOT::TGenericClassInfo 
         instance("unordered_map<string,TH2F*>", -2, "unordered_map", 1029,
                  typeid(unordered_map<string,TH2F*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &unordered_maplEstringcOTH2FmUgR_Dictionary, isa_proxy, 0,
                  sizeof(unordered_map<string,TH2F*>) );
      instance.SetNew(&new_unordered_maplEstringcOTH2FmUgR);
      instance.SetNewArray(&newArray_unordered_maplEstringcOTH2FmUgR);
      instance.SetDelete(&delete_unordered_maplEstringcOTH2FmUgR);
      instance.SetDeleteArray(&deleteArray_unordered_maplEstringcOTH2FmUgR);
      instance.SetDestructor(&destruct_unordered_maplEstringcOTH2FmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< unordered_map<string,TH2F*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("unordered_map<string,TH2F*>","std::__1::unordered_map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, TH2F*, std::__1::hash<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::equal_to<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, TH2F*>>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const unordered_map<string,TH2F*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *unordered_maplEstringcOTH2FmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const unordered_map<string,TH2F*>*>(nullptr))->GetClass();
      unordered_maplEstringcOTH2FmUgR_TClassManip(theClass);
   return theClass;
   }

   static void unordered_maplEstringcOTH2FmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_unordered_maplEstringcOTH2FmUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,TH2F*> : new unordered_map<string,TH2F*>;
   }
   static void *newArray_unordered_maplEstringcOTH2FmUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,TH2F*>[nElements] : new unordered_map<string,TH2F*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_unordered_maplEstringcOTH2FmUgR(void *p) {
      delete (static_cast<unordered_map<string,TH2F*>*>(p));
   }
   static void deleteArray_unordered_maplEstringcOTH2FmUgR(void *p) {
      delete [] (static_cast<unordered_map<string,TH2F*>*>(p));
   }
   static void destruct_unordered_maplEstringcOTH2FmUgR(void *p) {
      typedef unordered_map<string,TH2F*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class unordered_map<string,TH2F*>

namespace ROOT {
   static TClass *unordered_maplEstringcOTH1mUgR_Dictionary();
   static void unordered_maplEstringcOTH1mUgR_TClassManip(TClass*);
   static void *new_unordered_maplEstringcOTH1mUgR(void *p = nullptr);
   static void *newArray_unordered_maplEstringcOTH1mUgR(Long_t size, void *p);
   static void delete_unordered_maplEstringcOTH1mUgR(void *p);
   static void deleteArray_unordered_maplEstringcOTH1mUgR(void *p);
   static void destruct_unordered_maplEstringcOTH1mUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const unordered_map<string,TH1*>*)
   {
      unordered_map<string,TH1*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(unordered_map<string,TH1*>));
      static ::ROOT::TGenericClassInfo 
         instance("unordered_map<string,TH1*>", -2, "unordered_map", 1029,
                  typeid(unordered_map<string,TH1*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &unordered_maplEstringcOTH1mUgR_Dictionary, isa_proxy, 0,
                  sizeof(unordered_map<string,TH1*>) );
      instance.SetNew(&new_unordered_maplEstringcOTH1mUgR);
      instance.SetNewArray(&newArray_unordered_maplEstringcOTH1mUgR);
      instance.SetDelete(&delete_unordered_maplEstringcOTH1mUgR);
      instance.SetDeleteArray(&deleteArray_unordered_maplEstringcOTH1mUgR);
      instance.SetDestructor(&destruct_unordered_maplEstringcOTH1mUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< unordered_map<string,TH1*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("unordered_map<string,TH1*>","std::__1::unordered_map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, TH1*, std::__1::hash<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::equal_to<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, TH1*>>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const unordered_map<string,TH1*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *unordered_maplEstringcOTH1mUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const unordered_map<string,TH1*>*>(nullptr))->GetClass();
      unordered_maplEstringcOTH1mUgR_TClassManip(theClass);
   return theClass;
   }

   static void unordered_maplEstringcOTH1mUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_unordered_maplEstringcOTH1mUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,TH1*> : new unordered_map<string,TH1*>;
   }
   static void *newArray_unordered_maplEstringcOTH1mUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,TH1*>[nElements] : new unordered_map<string,TH1*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_unordered_maplEstringcOTH1mUgR(void *p) {
      delete (static_cast<unordered_map<string,TH1*>*>(p));
   }
   static void deleteArray_unordered_maplEstringcOTH1mUgR(void *p) {
      delete [] (static_cast<unordered_map<string,TH1*>*>(p));
   }
   static void destruct_unordered_maplEstringcOTH1mUgR(void *p) {
      typedef unordered_map<string,TH1*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class unordered_map<string,TH1*>

namespace ROOT {
   static TClass *unordered_maplEstringcOTF1mUgR_Dictionary();
   static void unordered_maplEstringcOTF1mUgR_TClassManip(TClass*);
   static void *new_unordered_maplEstringcOTF1mUgR(void *p = nullptr);
   static void *newArray_unordered_maplEstringcOTF1mUgR(Long_t size, void *p);
   static void delete_unordered_maplEstringcOTF1mUgR(void *p);
   static void deleteArray_unordered_maplEstringcOTF1mUgR(void *p);
   static void destruct_unordered_maplEstringcOTF1mUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const unordered_map<string,TF1*>*)
   {
      unordered_map<string,TF1*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(unordered_map<string,TF1*>));
      static ::ROOT::TGenericClassInfo 
         instance("unordered_map<string,TF1*>", -2, "unordered_map", 1029,
                  typeid(unordered_map<string,TF1*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &unordered_maplEstringcOTF1mUgR_Dictionary, isa_proxy, 0,
                  sizeof(unordered_map<string,TF1*>) );
      instance.SetNew(&new_unordered_maplEstringcOTF1mUgR);
      instance.SetNewArray(&newArray_unordered_maplEstringcOTF1mUgR);
      instance.SetDelete(&delete_unordered_maplEstringcOTF1mUgR);
      instance.SetDeleteArray(&deleteArray_unordered_maplEstringcOTF1mUgR);
      instance.SetDestructor(&destruct_unordered_maplEstringcOTF1mUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< unordered_map<string,TF1*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("unordered_map<string,TF1*>","std::__1::unordered_map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, TF1*, std::__1::hash<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::equal_to<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, TF1*>>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const unordered_map<string,TF1*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *unordered_maplEstringcOTF1mUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const unordered_map<string,TF1*>*>(nullptr))->GetClass();
      unordered_maplEstringcOTF1mUgR_TClassManip(theClass);
   return theClass;
   }

   static void unordered_maplEstringcOTF1mUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_unordered_maplEstringcOTF1mUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,TF1*> : new unordered_map<string,TF1*>;
   }
   static void *newArray_unordered_maplEstringcOTF1mUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,TF1*>[nElements] : new unordered_map<string,TF1*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_unordered_maplEstringcOTF1mUgR(void *p) {
      delete (static_cast<unordered_map<string,TF1*>*>(p));
   }
   static void deleteArray_unordered_maplEstringcOTF1mUgR(void *p) {
      delete [] (static_cast<unordered_map<string,TF1*>*>(p));
   }
   static void destruct_unordered_maplEstringcOTF1mUgR(void *p) {
      typedef unordered_map<string,TF1*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class unordered_map<string,TF1*>

namespace ROOT {
   static TClass *unordered_maplEstringcOPi0QAcLcLFitPairgR_Dictionary();
   static void unordered_maplEstringcOPi0QAcLcLFitPairgR_TClassManip(TClass*);
   static void *new_unordered_maplEstringcOPi0QAcLcLFitPairgR(void *p = nullptr);
   static void *newArray_unordered_maplEstringcOPi0QAcLcLFitPairgR(Long_t size, void *p);
   static void delete_unordered_maplEstringcOPi0QAcLcLFitPairgR(void *p);
   static void deleteArray_unordered_maplEstringcOPi0QAcLcLFitPairgR(void *p);
   static void destruct_unordered_maplEstringcOPi0QAcLcLFitPairgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const unordered_map<string,Pi0QA::FitPair>*)
   {
      unordered_map<string,Pi0QA::FitPair> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(unordered_map<string,Pi0QA::FitPair>));
      static ::ROOT::TGenericClassInfo 
         instance("unordered_map<string,Pi0QA::FitPair>", -2, "unordered_map", 1029,
                  typeid(unordered_map<string,Pi0QA::FitPair>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &unordered_maplEstringcOPi0QAcLcLFitPairgR_Dictionary, isa_proxy, 0,
                  sizeof(unordered_map<string,Pi0QA::FitPair>) );
      instance.SetNew(&new_unordered_maplEstringcOPi0QAcLcLFitPairgR);
      instance.SetNewArray(&newArray_unordered_maplEstringcOPi0QAcLcLFitPairgR);
      instance.SetDelete(&delete_unordered_maplEstringcOPi0QAcLcLFitPairgR);
      instance.SetDeleteArray(&deleteArray_unordered_maplEstringcOPi0QAcLcLFitPairgR);
      instance.SetDestructor(&destruct_unordered_maplEstringcOPi0QAcLcLFitPairgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< unordered_map<string,Pi0QA::FitPair> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("unordered_map<string,Pi0QA::FitPair>","std::__1::unordered_map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, Pi0QA::FitPair, std::__1::hash<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::equal_to<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, Pi0QA::FitPair>>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const unordered_map<string,Pi0QA::FitPair>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *unordered_maplEstringcOPi0QAcLcLFitPairgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const unordered_map<string,Pi0QA::FitPair>*>(nullptr))->GetClass();
      unordered_maplEstringcOPi0QAcLcLFitPairgR_TClassManip(theClass);
   return theClass;
   }

   static void unordered_maplEstringcOPi0QAcLcLFitPairgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_unordered_maplEstringcOPi0QAcLcLFitPairgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,Pi0QA::FitPair> : new unordered_map<string,Pi0QA::FitPair>;
   }
   static void *newArray_unordered_maplEstringcOPi0QAcLcLFitPairgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,Pi0QA::FitPair>[nElements] : new unordered_map<string,Pi0QA::FitPair>[nElements];
   }
   // Wrapper around operator delete
   static void delete_unordered_maplEstringcOPi0QAcLcLFitPairgR(void *p) {
      delete (static_cast<unordered_map<string,Pi0QA::FitPair>*>(p));
   }
   static void deleteArray_unordered_maplEstringcOPi0QAcLcLFitPairgR(void *p) {
      delete [] (static_cast<unordered_map<string,Pi0QA::FitPair>*>(p));
   }
   static void destruct_unordered_maplEstringcOPi0QAcLcLFitPairgR(void *p) {
      typedef unordered_map<string,Pi0QA::FitPair> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class unordered_map<string,Pi0QA::FitPair>

namespace ROOT {
   static TClass *unordered_maplEstringcOPi0QAcLcLFitInfogR_Dictionary();
   static void unordered_maplEstringcOPi0QAcLcLFitInfogR_TClassManip(TClass*);
   static void *new_unordered_maplEstringcOPi0QAcLcLFitInfogR(void *p = nullptr);
   static void *newArray_unordered_maplEstringcOPi0QAcLcLFitInfogR(Long_t size, void *p);
   static void delete_unordered_maplEstringcOPi0QAcLcLFitInfogR(void *p);
   static void deleteArray_unordered_maplEstringcOPi0QAcLcLFitInfogR(void *p);
   static void destruct_unordered_maplEstringcOPi0QAcLcLFitInfogR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const unordered_map<string,Pi0QA::FitInfo>*)
   {
      unordered_map<string,Pi0QA::FitInfo> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(unordered_map<string,Pi0QA::FitInfo>));
      static ::ROOT::TGenericClassInfo 
         instance("unordered_map<string,Pi0QA::FitInfo>", -2, "unordered_map", 1029,
                  typeid(unordered_map<string,Pi0QA::FitInfo>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &unordered_maplEstringcOPi0QAcLcLFitInfogR_Dictionary, isa_proxy, 0,
                  sizeof(unordered_map<string,Pi0QA::FitInfo>) );
      instance.SetNew(&new_unordered_maplEstringcOPi0QAcLcLFitInfogR);
      instance.SetNewArray(&newArray_unordered_maplEstringcOPi0QAcLcLFitInfogR);
      instance.SetDelete(&delete_unordered_maplEstringcOPi0QAcLcLFitInfogR);
      instance.SetDeleteArray(&deleteArray_unordered_maplEstringcOPi0QAcLcLFitInfogR);
      instance.SetDestructor(&destruct_unordered_maplEstringcOPi0QAcLcLFitInfogR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< unordered_map<string,Pi0QA::FitInfo> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("unordered_map<string,Pi0QA::FitInfo>","std::__1::unordered_map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, Pi0QA::FitInfo, std::__1::hash<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::equal_to<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, Pi0QA::FitInfo>>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const unordered_map<string,Pi0QA::FitInfo>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *unordered_maplEstringcOPi0QAcLcLFitInfogR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const unordered_map<string,Pi0QA::FitInfo>*>(nullptr))->GetClass();
      unordered_maplEstringcOPi0QAcLcLFitInfogR_TClassManip(theClass);
   return theClass;
   }

   static void unordered_maplEstringcOPi0QAcLcLFitInfogR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_unordered_maplEstringcOPi0QAcLcLFitInfogR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,Pi0QA::FitInfo> : new unordered_map<string,Pi0QA::FitInfo>;
   }
   static void *newArray_unordered_maplEstringcOPi0QAcLcLFitInfogR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,Pi0QA::FitInfo>[nElements] : new unordered_map<string,Pi0QA::FitInfo>[nElements];
   }
   // Wrapper around operator delete
   static void delete_unordered_maplEstringcOPi0QAcLcLFitInfogR(void *p) {
      delete (static_cast<unordered_map<string,Pi0QA::FitInfo>*>(p));
   }
   static void deleteArray_unordered_maplEstringcOPi0QAcLcLFitInfogR(void *p) {
      delete [] (static_cast<unordered_map<string,Pi0QA::FitInfo>*>(p));
   }
   static void destruct_unordered_maplEstringcOPi0QAcLcLFitInfogR(void *p) {
      typedef unordered_map<string,Pi0QA::FitInfo> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class unordered_map<string,Pi0QA::FitInfo>

namespace ROOT {
   static TClass *unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR_Dictionary();
   static void unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR_TClassManip(TClass*);
   static void *new_unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR(void *p = nullptr);
   static void *newArray_unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR(Long_t size, void *p);
   static void delete_unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR(void *p);
   static void deleteArray_unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR(void *p);
   static void destruct_unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const unordered_map<string,MapPair,NSCache<MapPair>::H>*)
   {
      unordered_map<string,MapPair,NSCache<MapPair>::H> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(unordered_map<string,MapPair,NSCache<MapPair>::H>));
      static ::ROOT::TGenericClassInfo 
         instance("unordered_map<string,MapPair,NSCache<MapPair>::H>", -2, "unordered_map", 1029,
                  typeid(unordered_map<string,MapPair,NSCache<MapPair>::H>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR_Dictionary, isa_proxy, 0,
                  sizeof(unordered_map<string,MapPair,NSCache<MapPair>::H>) );
      instance.SetNew(&new_unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR);
      instance.SetNewArray(&newArray_unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR);
      instance.SetDelete(&delete_unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR);
      instance.SetDeleteArray(&deleteArray_unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR);
      instance.SetDestructor(&destruct_unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< unordered_map<string,MapPair,NSCache<MapPair>::H> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("unordered_map<string,MapPair,NSCache<MapPair>::H>","std::__1::unordered_map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, MapPair, NSCache<MapPair>::H, std::__1::equal_to<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, MapPair>>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const unordered_map<string,MapPair,NSCache<MapPair>::H>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const unordered_map<string,MapPair,NSCache<MapPair>::H>*>(nullptr))->GetClass();
      unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR_TClassManip(theClass);
   return theClass;
   }

   static void unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,MapPair,NSCache<MapPair>::H> : new unordered_map<string,MapPair,NSCache<MapPair>::H>;
   }
   static void *newArray_unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) unordered_map<string,MapPair,NSCache<MapPair>::H>[nElements] : new unordered_map<string,MapPair,NSCache<MapPair>::H>[nElements];
   }
   // Wrapper around operator delete
   static void delete_unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR(void *p) {
      delete (static_cast<unordered_map<string,MapPair,NSCache<MapPair>::H>*>(p));
   }
   static void deleteArray_unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR(void *p) {
      delete [] (static_cast<unordered_map<string,MapPair,NSCache<MapPair>::H>*>(p));
   }
   static void destruct_unordered_maplEstringcOMapPaircONSCachelEMapPairgRcLcLHgR(void *p) {
      typedef unordered_map<string,MapPair,NSCache<MapPair>::H> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class unordered_map<string,MapPair,NSCache<MapPair>::H>

namespace ROOT {
   static TClass *maplEstringcOvectorlETProfilemUgRsPgR_Dictionary();
   static void maplEstringcOvectorlETProfilemUgRsPgR_TClassManip(TClass*);
   static void *new_maplEstringcOvectorlETProfilemUgRsPgR(void *p = nullptr);
   static void *newArray_maplEstringcOvectorlETProfilemUgRsPgR(Long_t size, void *p);
   static void delete_maplEstringcOvectorlETProfilemUgRsPgR(void *p);
   static void deleteArray_maplEstringcOvectorlETProfilemUgRsPgR(void *p);
   static void destruct_maplEstringcOvectorlETProfilemUgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,vector<TProfile*> >*)
   {
      map<string,vector<TProfile*> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,vector<TProfile*> >));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,vector<TProfile*> >", -2, "map", 967,
                  typeid(map<string,vector<TProfile*> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOvectorlETProfilemUgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,vector<TProfile*> >) );
      instance.SetNew(&new_maplEstringcOvectorlETProfilemUgRsPgR);
      instance.SetNewArray(&newArray_maplEstringcOvectorlETProfilemUgRsPgR);
      instance.SetDelete(&delete_maplEstringcOvectorlETProfilemUgRsPgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOvectorlETProfilemUgRsPgR);
      instance.SetDestructor(&destruct_maplEstringcOvectorlETProfilemUgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,vector<TProfile*> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,vector<TProfile*> >","std::__1::map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>, std::__1::less<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>>>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,vector<TProfile*> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOvectorlETProfilemUgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,vector<TProfile*> >*>(nullptr))->GetClass();
      maplEstringcOvectorlETProfilemUgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOvectorlETProfilemUgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOvectorlETProfilemUgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,vector<TProfile*> > : new map<string,vector<TProfile*> >;
   }
   static void *newArray_maplEstringcOvectorlETProfilemUgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,vector<TProfile*> >[nElements] : new map<string,vector<TProfile*> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOvectorlETProfilemUgRsPgR(void *p) {
      delete (static_cast<map<string,vector<TProfile*> >*>(p));
   }
   static void deleteArray_maplEstringcOvectorlETProfilemUgRsPgR(void *p) {
      delete [] (static_cast<map<string,vector<TProfile*> >*>(p));
   }
   static void destruct_maplEstringcOvectorlETProfilemUgRsPgR(void *p) {
      typedef map<string,vector<TProfile*> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,vector<TProfile*> >

namespace ROOT {
   static TClass *maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR_Dictionary();
   static void maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR_TClassManip(TClass*);
   static void *new_maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR(void *p = nullptr);
   static void *newArray_maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR(Long_t size, void *p);
   static void delete_maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR(void *p);
   static void deleteArray_maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR(void *p);
   static void destruct_maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,map<string,vector<TProfile*> > >*)
   {
      map<int,map<string,vector<TProfile*> > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,map<string,vector<TProfile*> > >));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,map<string,vector<TProfile*> > >", -2, "map", 967,
                  typeid(map<int,map<string,vector<TProfile*> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,map<string,vector<TProfile*> > >) );
      instance.SetNew(&new_maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR);
      instance.SetNewArray(&newArray_maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR);
      instance.SetDelete(&delete_maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR);
      instance.SetDestructor(&destruct_maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,map<string,vector<TProfile*> > > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<int,map<string,vector<TProfile*> > >","std::__1::map<int, std::__1::map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>, std::__1::less<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>>>>, std::__1::less<int>, std::__1::allocator<std::__1::pair<int const, std::__1::map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>, std::__1::less<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>>, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> const, std::__1::vector<TProfile*, std::__1::allocator<TProfile*>>>>>>>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<int,map<string,vector<TProfile*> > >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<int,map<string,vector<TProfile*> > >*>(nullptr))->GetClass();
      maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<int,map<string,vector<TProfile*> > > : new map<int,map<string,vector<TProfile*> > >;
   }
   static void *newArray_maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<int,map<string,vector<TProfile*> > >[nElements] : new map<int,map<string,vector<TProfile*> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR(void *p) {
      delete (static_cast<map<int,map<string,vector<TProfile*> > >*>(p));
   }
   static void deleteArray_maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR(void *p) {
      delete [] (static_cast<map<int,map<string,vector<TProfile*> > >*>(p));
   }
   static void destruct_maplEintcOmaplEstringcOvectorlETProfilemUgRsPgRsPgR(void *p) {
      typedef map<int,map<string,vector<TProfile*> > > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<int,map<string,vector<TProfile*> > >

namespace {
  void TriggerDictionaryInitialization_analyzeRun24or25auau_cpp_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./analyzeRun24or25auau.cpp",
nullptr
    };
    static const char* includePaths[] = {
"/Users/patsfan753/Desktop/analysis/env/include",
"/Users/patsfan753/Desktop/analysis/env/etc/",
"/Users/patsfan753/Desktop/analysis/env/etc//cling",
"/Users/patsfan753/Desktop/analysis/env/etc//cling/plugins/include",
"/Users/patsfan753/Desktop/analysis/env/include/",
"/Users/patsfan753/Desktop/analysis/env/include",
"/Users/runner/miniforge3/conda-bld/bld/rattler-build_root_base_1740585318/work/build-dir/include",
"/Users/runner/miniforge3/conda-bld/bld/rattler-build_root_base_1740585318/work/build-rootcling-xp/include/",
"/Users/patsfan753/Desktop/analysis/env/include/",
"/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/macros/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "analyzeRun24or25auau_cpp_ACLiC_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
struct __attribute__((annotate("$clingAutoload$./analyzeRun24or25auau.cpp")))  CutKey;
struct __attribute__((annotate("$clingAutoload$./analyzeRun24or25auau.cpp")))  MapPair;
template <class MPair> class __attribute__((annotate("$clingAutoload$./analyzeRun24or25auau.cpp")))  NSCache;

class __attribute__((annotate("$clingAutoload$./analyzeRun24or25auau.cpp")))  QA;
class __attribute__((annotate("$clingAutoload$./analyzeRun24or25auau.cpp")))  Pi0QA;
struct __attribute__((annotate("$clingAutoload$./analyzeRun24or25auau.cpp")))  NSPair;
class __attribute__((annotate("$clingAutoload$./analyzeRun24or25auau.cpp")))  CorrQA;
class __attribute__((annotate("$clingAutoload$./analyzeRun24or25auau.cpp")))  EmcalQA;
class __attribute__((annotate("$clingAutoload$./analyzeRun24or25auau.cpp")))  HcalQA;
class __attribute__((annotate("$clingAutoload$./analyzeRun24or25auau.cpp")))  SepdPlaneQA;
struct __attribute__((annotate("$clingAutoload$./analyzeRun24or25auau.cpp")))  MBDTag;
template <class DERIVED> class __attribute__((annotate("$clingAutoload$./analyzeRun24or25auau.cpp")))  NSDetectorQA;

struct __attribute__((annotate("$clingAutoload$./analyzeRun24or25auau.cpp")))  sEPDTag;
class __attribute__((annotate("$clingAutoload$./analyzeRun24or25auau.cpp")))  EventQA;
class __attribute__((annotate("$clingAutoload$./analyzeRun24or25auau.cpp")))  JetQA;
class __attribute__((annotate("$clingAutoload$./analyzeRun24or25auau.cpp")))  VnPlotQA;
using MbdQA __attribute__((annotate("$clingAutoload$./analyzeRun24or25auau.cpp")))  = NSDetectorQA<MBDTag>;
using SepdQA __attribute__((annotate("$clingAutoload$./analyzeRun24or25auau.cpp")))  = NSDetectorQA<sEPDTag>;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "analyzeRun24or25auau_cpp_ACLiC_dict dictionary payload"

#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "./analyzeRun24or25auau.cpp"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"", payloadCode, "@",
"CentList", payloadCode, "@",
"CorrQA", payloadCode, "@",
"CorrQA::s_cache", payloadCode, "@",
"CutKey", payloadCode, "@",
"EmcalQA", payloadCode, "@",
"EmcalQA::kFracSaturation", payloadCode, "@",
"EventQA", payloadCode, "@",
"EventQA::s_centHists", payloadCode, "@",
"EventQA::s_evtCounts", payloadCode, "@",
"EventQA::s_points", payloadCode, "@",
"EventQA::s_summaryWritten", payloadCode, "@",
"HcalQA", payloadCode, "@",
"JetQA", payloadCode, "@",
"MBDTag", payloadCode, "@",
"MBDTag::subdir", payloadCode, "@",
"MBDTag::titleNorth", payloadCode, "@",
"MBDTag::titleSouth", payloadCode, "@",
"MapPair", payloadCode, "@",
"MbdQA", payloadCode, "@",
"NSCache<MapPair>", payloadCode, "@",
"NSDetectorQA<MBDTag>", payloadCode, "@",
"NSDetectorQA<sEPDTag>", payloadCode, "@",
"NSPair", payloadCode, "@",
"Pi0QA", payloadCode, "@",
"Pi0QA::s_runPoints", payloadCode, "@",
"Pi0QA::s_summaryWritten", payloadCode, "@",
"QA", payloadCode, "@",
"SepdPlaneQA", payloadCode, "@",
"SepdQA", payloadCode, "@",
"VnPlotQA", payloadCode, "@",
"analyzeRun24or25auau", payloadCode, "@",
"binAt", payloadCode, "@",
"cPath", payloadCode, "@",
"decodeInvName", payloadCode, "@",
"discoverSlices", payloadCode, "@",
"drawRunLabel", payloadCode, "@",
"ensure_dir", payloadCode, "@",
"g_nsCache", payloadCode, "@",
"hcal_plate_from_idx", payloadCode, "@",
"hcal_sector_from_idx", payloadCode, "@",
"ib_from_idx", payloadCode, "@",
"isBadBoard", payloadCode, "@",
"isBadHcalPlate", payloadCode, "@",
"kDoPi0Fit", payloadCode, "@",
"kEMCalCanvasH", payloadCode, "@",
"kEMCalCanvasW", payloadCode, "@",
"kHCalCanvasH", payloadCode, "@",
"kHCalCanvasW", payloadCode, "@",
"kInputFile", payloadCode, "@",
"kOutputBase", payloadCode, "@",
"kTriggersWanted", payloadCode, "@",
"listRunFiles", payloadCode, "@",
"log::banner", payloadCode, "@",
"log::err", payloadCode, "@",
"log::info", payloadCode, "@",
"log::ok", payloadCode, "@",
"log::trace", payloadCode, "@",
"log::warn", payloadCode, "@",
"runOneQaPass", payloadCode, "@",
"sEPDTag", payloadCode, "@",
"sEPDTag::subdir", payloadCode, "@",
"sEPDTag::titleNorth", payloadCode, "@",
"sEPDTag::titleSouth", payloadCode, "@",
"save1D", payloadCode, "@",
"save2D", payloadCode, "@",
"sector_from_idx", payloadCode, "@",
"setupPad", payloadCode, "@",
"sf3", payloadCode, "@",
"stripLeadingZeros", payloadCode, "@",
"styleAxes", payloadCode, "@",
"term::CLR_BOLD", payloadCode, "@",
"term::CLR_CYAN", payloadCode, "@",
"term::CLR_GRN", payloadCode, "@",
"term::CLR_RED", payloadCode, "@",
"term::CLR_RST", payloadCode, "@",
"term::CLR_YEL", payloadCode, "@",
"tidyAxes", payloadCode, "@",
"tightenAxes", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("analyzeRun24or25auau_cpp_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_analyzeRun24or25auau_cpp_ACLiC_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_analyzeRun24or25auau_cpp_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_analyzeRun24or25auau_cpp_ACLiC_dict() {
  TriggerDictionaryInitialization_analyzeRun24or25auau_cpp_ACLiC_dict_Impl();
}
