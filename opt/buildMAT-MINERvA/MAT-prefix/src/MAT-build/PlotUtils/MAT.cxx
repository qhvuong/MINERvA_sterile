// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME MAT
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

// Header files passed as explicit arguments
#include "/exp/minerva/app/users/qvuong/MAT_AL9/MAT/PlotUtils/../dict/PlotUtilsDict.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *PlotUtilscLcLBadUnivWgtError_Dictionary();
   static void PlotUtilscLcLBadUnivWgtError_TClassManip(TClass*);
   static void delete_PlotUtilscLcLBadUnivWgtError(void *p);
   static void deleteArray_PlotUtilscLcLBadUnivWgtError(void *p);
   static void destruct_PlotUtilscLcLBadUnivWgtError(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::BadUnivWgtError*)
   {
      ::PlotUtils::BadUnivWgtError *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::PlotUtils::BadUnivWgtError));
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::BadUnivWgtError", "Exceptions.h", 16,
                  typeid(::PlotUtils::BadUnivWgtError), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &PlotUtilscLcLBadUnivWgtError_Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::BadUnivWgtError) );
      instance.SetDelete(&delete_PlotUtilscLcLBadUnivWgtError);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLBadUnivWgtError);
      instance.SetDestructor(&destruct_PlotUtilscLcLBadUnivWgtError);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::BadUnivWgtError*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::BadUnivWgtError*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::BadUnivWgtError*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *PlotUtilscLcLBadUnivWgtError_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::PlotUtils::BadUnivWgtError*>(nullptr))->GetClass();
      PlotUtilscLcLBadUnivWgtError_TClassManip(theClass);
   return theClass;
   }

   static void PlotUtilscLcLBadUnivWgtError_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *PlotUtilscLcLConstraintLoadError_Dictionary();
   static void PlotUtilscLcLConstraintLoadError_TClassManip(TClass*);
   static void delete_PlotUtilscLcLConstraintLoadError(void *p);
   static void deleteArray_PlotUtilscLcLConstraintLoadError(void *p);
   static void destruct_PlotUtilscLcLConstraintLoadError(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::ConstraintLoadError*)
   {
      ::PlotUtils::ConstraintLoadError *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::PlotUtils::ConstraintLoadError));
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::ConstraintLoadError", "Exceptions.h", 24,
                  typeid(::PlotUtils::ConstraintLoadError), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &PlotUtilscLcLConstraintLoadError_Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::ConstraintLoadError) );
      instance.SetDelete(&delete_PlotUtilscLcLConstraintLoadError);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLConstraintLoadError);
      instance.SetDestructor(&destruct_PlotUtilscLcLConstraintLoadError);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::ConstraintLoadError*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::ConstraintLoadError*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::ConstraintLoadError*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *PlotUtilscLcLConstraintLoadError_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::PlotUtils::ConstraintLoadError*>(nullptr))->GetClass();
      PlotUtilscLcLConstraintLoadError_TClassManip(theClass);
   return theClass;
   }

   static void PlotUtilscLcLConstraintLoadError_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *PlotUtilscLcLConstraintAccessError_Dictionary();
   static void PlotUtilscLcLConstraintAccessError_TClassManip(TClass*);
   static void delete_PlotUtilscLcLConstraintAccessError(void *p);
   static void deleteArray_PlotUtilscLcLConstraintAccessError(void *p);
   static void destruct_PlotUtilscLcLConstraintAccessError(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::ConstraintAccessError*)
   {
      ::PlotUtils::ConstraintAccessError *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::PlotUtils::ConstraintAccessError));
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::ConstraintAccessError", "Exceptions.h", 32,
                  typeid(::PlotUtils::ConstraintAccessError), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &PlotUtilscLcLConstraintAccessError_Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::ConstraintAccessError) );
      instance.SetDelete(&delete_PlotUtilscLcLConstraintAccessError);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLConstraintAccessError);
      instance.SetDestructor(&destruct_PlotUtilscLcLConstraintAccessError);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::ConstraintAccessError*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::ConstraintAccessError*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::ConstraintAccessError*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *PlotUtilscLcLConstraintAccessError_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::PlotUtils::ConstraintAccessError*>(nullptr))->GetClass();
      PlotUtilscLcLConstraintAccessError_TClassManip(theClass);
   return theClass;
   }

   static void PlotUtilscLcLConstraintAccessError_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *PlotUtilscLcLMissingFluxUnivError_Dictionary();
   static void PlotUtilscLcLMissingFluxUnivError_TClassManip(TClass*);
   static void delete_PlotUtilscLcLMissingFluxUnivError(void *p);
   static void deleteArray_PlotUtilscLcLMissingFluxUnivError(void *p);
   static void destruct_PlotUtilscLcLMissingFluxUnivError(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::MissingFluxUnivError*)
   {
      ::PlotUtils::MissingFluxUnivError *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::PlotUtils::MissingFluxUnivError));
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::MissingFluxUnivError", "Exceptions.h", 40,
                  typeid(::PlotUtils::MissingFluxUnivError), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &PlotUtilscLcLMissingFluxUnivError_Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::MissingFluxUnivError) );
      instance.SetDelete(&delete_PlotUtilscLcLMissingFluxUnivError);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLMissingFluxUnivError);
      instance.SetDestructor(&destruct_PlotUtilscLcLMissingFluxUnivError);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::MissingFluxUnivError*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::MissingFluxUnivError*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::MissingFluxUnivError*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *PlotUtilscLcLMissingFluxUnivError_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::PlotUtils::MissingFluxUnivError*>(nullptr))->GetClass();
      PlotUtilscLcLMissingFluxUnivError_TClassManip(theClass);
   return theClass;
   }

   static void PlotUtilscLcLMissingFluxUnivError_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *PlotUtilscLcLMissingSpectatorStrategyError_Dictionary();
   static void PlotUtilscLcLMissingSpectatorStrategyError_TClassManip(TClass*);
   static void delete_PlotUtilscLcLMissingSpectatorStrategyError(void *p);
   static void deleteArray_PlotUtilscLcLMissingSpectatorStrategyError(void *p);
   static void destruct_PlotUtilscLcLMissingSpectatorStrategyError(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::MissingSpectatorStrategyError*)
   {
      ::PlotUtils::MissingSpectatorStrategyError *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::PlotUtils::MissingSpectatorStrategyError));
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::MissingSpectatorStrategyError", "Exceptions.h", 48,
                  typeid(::PlotUtils::MissingSpectatorStrategyError), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &PlotUtilscLcLMissingSpectatorStrategyError_Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::MissingSpectatorStrategyError) );
      instance.SetDelete(&delete_PlotUtilscLcLMissingSpectatorStrategyError);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLMissingSpectatorStrategyError);
      instance.SetDestructor(&destruct_PlotUtilscLcLMissingSpectatorStrategyError);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::MissingSpectatorStrategyError*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::MissingSpectatorStrategyError*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::MissingSpectatorStrategyError*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *PlotUtilscLcLMissingSpectatorStrategyError_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::PlotUtils::MissingSpectatorStrategyError*>(nullptr))->GetClass();
      PlotUtilscLcLMissingSpectatorStrategyError_TClassManip(theClass);
   return theClass;
   }

   static void PlotUtilscLcLMissingSpectatorStrategyError_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *PlotUtilscLcLNoFluxUnivError_Dictionary();
   static void PlotUtilscLcLNoFluxUnivError_TClassManip(TClass*);
   static void delete_PlotUtilscLcLNoFluxUnivError(void *p);
   static void deleteArray_PlotUtilscLcLNoFluxUnivError(void *p);
   static void destruct_PlotUtilscLcLNoFluxUnivError(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::NoFluxUnivError*)
   {
      ::PlotUtils::NoFluxUnivError *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::PlotUtils::NoFluxUnivError));
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::NoFluxUnivError", "Exceptions.h", 56,
                  typeid(::PlotUtils::NoFluxUnivError), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &PlotUtilscLcLNoFluxUnivError_Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::NoFluxUnivError) );
      instance.SetDelete(&delete_PlotUtilscLcLNoFluxUnivError);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLNoFluxUnivError);
      instance.SetDestructor(&destruct_PlotUtilscLcLNoFluxUnivError);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::NoFluxUnivError*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::NoFluxUnivError*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::NoFluxUnivError*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *PlotUtilscLcLNoFluxUnivError_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::PlotUtils::NoFluxUnivError*>(nullptr))->GetClass();
      PlotUtilscLcLNoFluxUnivError_TClassManip(theClass);
   return theClass;
   }

   static void PlotUtilscLcLNoFluxUnivError_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *PlotUtilscLcLNoWgtdSpreadError_Dictionary();
   static void PlotUtilscLcLNoWgtdSpreadError_TClassManip(TClass*);
   static void delete_PlotUtilscLcLNoWgtdSpreadError(void *p);
   static void deleteArray_PlotUtilscLcLNoWgtdSpreadError(void *p);
   static void destruct_PlotUtilscLcLNoWgtdSpreadError(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::NoWgtdSpreadError*)
   {
      ::PlotUtils::NoWgtdSpreadError *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::PlotUtils::NoWgtdSpreadError));
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::NoWgtdSpreadError", "Exceptions.h", 65,
                  typeid(::PlotUtils::NoWgtdSpreadError), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &PlotUtilscLcLNoWgtdSpreadError_Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::NoWgtdSpreadError) );
      instance.SetDelete(&delete_PlotUtilscLcLNoWgtdSpreadError);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLNoWgtdSpreadError);
      instance.SetDestructor(&destruct_PlotUtilscLcLNoWgtdSpreadError);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::NoWgtdSpreadError*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::NoWgtdSpreadError*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::NoWgtdSpreadError*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *PlotUtilscLcLNoWgtdSpreadError_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::PlotUtils::NoWgtdSpreadError*>(nullptr))->GetClass();
      PlotUtilscLcLNoWgtdSpreadError_TClassManip(theClass);
   return theClass;
   }

   static void PlotUtilscLcLNoWgtdSpreadError_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *PlotUtilscLcLSpectatorConstraintCollisionError_Dictionary();
   static void PlotUtilscLcLSpectatorConstraintCollisionError_TClassManip(TClass*);
   static void delete_PlotUtilscLcLSpectatorConstraintCollisionError(void *p);
   static void deleteArray_PlotUtilscLcLSpectatorConstraintCollisionError(void *p);
   static void destruct_PlotUtilscLcLSpectatorConstraintCollisionError(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::SpectatorConstraintCollisionError*)
   {
      ::PlotUtils::SpectatorConstraintCollisionError *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::PlotUtils::SpectatorConstraintCollisionError));
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::SpectatorConstraintCollisionError", "Exceptions.h", 73,
                  typeid(::PlotUtils::SpectatorConstraintCollisionError), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &PlotUtilscLcLSpectatorConstraintCollisionError_Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::SpectatorConstraintCollisionError) );
      instance.SetDelete(&delete_PlotUtilscLcLSpectatorConstraintCollisionError);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLSpectatorConstraintCollisionError);
      instance.SetDestructor(&destruct_PlotUtilscLcLSpectatorConstraintCollisionError);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::SpectatorConstraintCollisionError*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::SpectatorConstraintCollisionError*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::SpectatorConstraintCollisionError*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *PlotUtilscLcLSpectatorConstraintCollisionError_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::PlotUtils::SpectatorConstraintCollisionError*>(nullptr))->GetClass();
      PlotUtilscLcLSpectatorConstraintCollisionError_TClassManip(theClass);
   return theClass;
   }

   static void PlotUtilscLcLSpectatorConstraintCollisionError_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_PlotUtilscLcLMnvVertErrorBand(void *p = nullptr);
   static void *newArray_PlotUtilscLcLMnvVertErrorBand(Long_t size, void *p);
   static void delete_PlotUtilscLcLMnvVertErrorBand(void *p);
   static void deleteArray_PlotUtilscLcLMnvVertErrorBand(void *p);
   static void destruct_PlotUtilscLcLMnvVertErrorBand(void *p);
   static void directoryAutoAdd_PlotUtilscLcLMnvVertErrorBand(void *obj, TDirectory *dir);
   static Long64_t merge_PlotUtilscLcLMnvVertErrorBand(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::MnvVertErrorBand*)
   {
      ::PlotUtils::MnvVertErrorBand *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PlotUtils::MnvVertErrorBand >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::MnvVertErrorBand", ::PlotUtils::MnvVertErrorBand::Class_Version(), "MnvVertErrorBand.h", 33,
                  typeid(::PlotUtils::MnvVertErrorBand), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PlotUtils::MnvVertErrorBand::Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::MnvVertErrorBand) );
      instance.SetNew(&new_PlotUtilscLcLMnvVertErrorBand);
      instance.SetNewArray(&newArray_PlotUtilscLcLMnvVertErrorBand);
      instance.SetDelete(&delete_PlotUtilscLcLMnvVertErrorBand);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLMnvVertErrorBand);
      instance.SetDestructor(&destruct_PlotUtilscLcLMnvVertErrorBand);
      instance.SetDirectoryAutoAdd(&directoryAutoAdd_PlotUtilscLcLMnvVertErrorBand);
      instance.SetMerge(&merge_PlotUtilscLcLMnvVertErrorBand);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::MnvVertErrorBand*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::MnvVertErrorBand*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::MnvVertErrorBand*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PlotUtilscLcLMnvLatErrorBand(void *p = nullptr);
   static void *newArray_PlotUtilscLcLMnvLatErrorBand(Long_t size, void *p);
   static void delete_PlotUtilscLcLMnvLatErrorBand(void *p);
   static void deleteArray_PlotUtilscLcLMnvLatErrorBand(void *p);
   static void destruct_PlotUtilscLcLMnvLatErrorBand(void *p);
   static void directoryAutoAdd_PlotUtilscLcLMnvLatErrorBand(void *obj, TDirectory *dir);
   static Long64_t merge_PlotUtilscLcLMnvLatErrorBand(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::MnvLatErrorBand*)
   {
      ::PlotUtils::MnvLatErrorBand *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PlotUtils::MnvLatErrorBand >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::MnvLatErrorBand", ::PlotUtils::MnvLatErrorBand::Class_Version(), "MnvLatErrorBand.h", 34,
                  typeid(::PlotUtils::MnvLatErrorBand), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PlotUtils::MnvLatErrorBand::Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::MnvLatErrorBand) );
      instance.SetNew(&new_PlotUtilscLcLMnvLatErrorBand);
      instance.SetNewArray(&newArray_PlotUtilscLcLMnvLatErrorBand);
      instance.SetDelete(&delete_PlotUtilscLcLMnvLatErrorBand);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLMnvLatErrorBand);
      instance.SetDestructor(&destruct_PlotUtilscLcLMnvLatErrorBand);
      instance.SetDirectoryAutoAdd(&directoryAutoAdd_PlotUtilscLcLMnvLatErrorBand);
      instance.SetMerge(&merge_PlotUtilscLcLMnvLatErrorBand);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::MnvLatErrorBand*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::MnvLatErrorBand*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::MnvLatErrorBand*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PlotUtilscLcLMnvH1D(void *p = nullptr);
   static void *newArray_PlotUtilscLcLMnvH1D(Long_t size, void *p);
   static void delete_PlotUtilscLcLMnvH1D(void *p);
   static void deleteArray_PlotUtilscLcLMnvH1D(void *p);
   static void destruct_PlotUtilscLcLMnvH1D(void *p);
   static void directoryAutoAdd_PlotUtilscLcLMnvH1D(void *obj, TDirectory *dir);
   static Long64_t merge_PlotUtilscLcLMnvH1D(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::MnvH1D*)
   {
      ::PlotUtils::MnvH1D *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PlotUtils::MnvH1D >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::MnvH1D", ::PlotUtils::MnvH1D::Class_Version(), "MnvH1D.h", 41,
                  typeid(::PlotUtils::MnvH1D), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PlotUtils::MnvH1D::Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::MnvH1D) );
      instance.SetNew(&new_PlotUtilscLcLMnvH1D);
      instance.SetNewArray(&newArray_PlotUtilscLcLMnvH1D);
      instance.SetDelete(&delete_PlotUtilscLcLMnvH1D);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLMnvH1D);
      instance.SetDestructor(&destruct_PlotUtilscLcLMnvH1D);
      instance.SetDirectoryAutoAdd(&directoryAutoAdd_PlotUtilscLcLMnvH1D);
      instance.SetMerge(&merge_PlotUtilscLcLMnvH1D);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::MnvH1D*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::MnvH1D*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::MnvH1D*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *PlotUtilscLcLMnvHist_Dictionary();
   static void PlotUtilscLcLMnvHist_TClassManip(TClass*);
   static void *new_PlotUtilscLcLMnvHist(void *p = nullptr);
   static void *newArray_PlotUtilscLcLMnvHist(Long_t size, void *p);
   static void delete_PlotUtilscLcLMnvHist(void *p);
   static void deleteArray_PlotUtilscLcLMnvHist(void *p);
   static void destruct_PlotUtilscLcLMnvHist(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::MnvHist*)
   {
      ::PlotUtils::MnvHist *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::PlotUtils::MnvHist));
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::MnvHist", "HistogramUtils.h", 33,
                  typeid(::PlotUtils::MnvHist), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &PlotUtilscLcLMnvHist_Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::MnvHist) );
      instance.SetNew(&new_PlotUtilscLcLMnvHist);
      instance.SetNewArray(&newArray_PlotUtilscLcLMnvHist);
      instance.SetDelete(&delete_PlotUtilscLcLMnvHist);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLMnvHist);
      instance.SetDestructor(&destruct_PlotUtilscLcLMnvHist);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::MnvHist*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::MnvHist*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::MnvHist*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *PlotUtilscLcLMnvHist_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::PlotUtils::MnvHist*>(nullptr))->GetClass();
      PlotUtilscLcLMnvHist_TClassManip(theClass);
   return theClass;
   }

   static void PlotUtilscLcLMnvHist_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_PlotUtilscLcLMnvVertErrorBand2D(void *p = nullptr);
   static void *newArray_PlotUtilscLcLMnvVertErrorBand2D(Long_t size, void *p);
   static void delete_PlotUtilscLcLMnvVertErrorBand2D(void *p);
   static void deleteArray_PlotUtilscLcLMnvVertErrorBand2D(void *p);
   static void destruct_PlotUtilscLcLMnvVertErrorBand2D(void *p);
   static void directoryAutoAdd_PlotUtilscLcLMnvVertErrorBand2D(void *obj, TDirectory *dir);
   static Long64_t merge_PlotUtilscLcLMnvVertErrorBand2D(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::MnvVertErrorBand2D*)
   {
      ::PlotUtils::MnvVertErrorBand2D *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PlotUtils::MnvVertErrorBand2D >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::MnvVertErrorBand2D", ::PlotUtils::MnvVertErrorBand2D::Class_Version(), "MnvVertErrorBand2D.h", 21,
                  typeid(::PlotUtils::MnvVertErrorBand2D), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PlotUtils::MnvVertErrorBand2D::Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::MnvVertErrorBand2D) );
      instance.SetNew(&new_PlotUtilscLcLMnvVertErrorBand2D);
      instance.SetNewArray(&newArray_PlotUtilscLcLMnvVertErrorBand2D);
      instance.SetDelete(&delete_PlotUtilscLcLMnvVertErrorBand2D);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLMnvVertErrorBand2D);
      instance.SetDestructor(&destruct_PlotUtilscLcLMnvVertErrorBand2D);
      instance.SetDirectoryAutoAdd(&directoryAutoAdd_PlotUtilscLcLMnvVertErrorBand2D);
      instance.SetMerge(&merge_PlotUtilscLcLMnvVertErrorBand2D);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::MnvVertErrorBand2D*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::MnvVertErrorBand2D*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::MnvVertErrorBand2D*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PlotUtilscLcLMnvLatErrorBand2D(void *p = nullptr);
   static void *newArray_PlotUtilscLcLMnvLatErrorBand2D(Long_t size, void *p);
   static void delete_PlotUtilscLcLMnvLatErrorBand2D(void *p);
   static void deleteArray_PlotUtilscLcLMnvLatErrorBand2D(void *p);
   static void destruct_PlotUtilscLcLMnvLatErrorBand2D(void *p);
   static void directoryAutoAdd_PlotUtilscLcLMnvLatErrorBand2D(void *obj, TDirectory *dir);
   static Long64_t merge_PlotUtilscLcLMnvLatErrorBand2D(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::MnvLatErrorBand2D*)
   {
      ::PlotUtils::MnvLatErrorBand2D *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PlotUtils::MnvLatErrorBand2D >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::MnvLatErrorBand2D", ::PlotUtils::MnvLatErrorBand2D::Class_Version(), "MnvLatErrorBand2D.h", 21,
                  typeid(::PlotUtils::MnvLatErrorBand2D), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PlotUtils::MnvLatErrorBand2D::Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::MnvLatErrorBand2D) );
      instance.SetNew(&new_PlotUtilscLcLMnvLatErrorBand2D);
      instance.SetNewArray(&newArray_PlotUtilscLcLMnvLatErrorBand2D);
      instance.SetDelete(&delete_PlotUtilscLcLMnvLatErrorBand2D);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLMnvLatErrorBand2D);
      instance.SetDestructor(&destruct_PlotUtilscLcLMnvLatErrorBand2D);
      instance.SetDirectoryAutoAdd(&directoryAutoAdd_PlotUtilscLcLMnvLatErrorBand2D);
      instance.SetMerge(&merge_PlotUtilscLcLMnvLatErrorBand2D);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::MnvLatErrorBand2D*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::MnvLatErrorBand2D*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::MnvLatErrorBand2D*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PlotUtilscLcLMnvH2D(void *p = nullptr);
   static void *newArray_PlotUtilscLcLMnvH2D(Long_t size, void *p);
   static void delete_PlotUtilscLcLMnvH2D(void *p);
   static void deleteArray_PlotUtilscLcLMnvH2D(void *p);
   static void destruct_PlotUtilscLcLMnvH2D(void *p);
   static void directoryAutoAdd_PlotUtilscLcLMnvH2D(void *obj, TDirectory *dir);
   static Long64_t merge_PlotUtilscLcLMnvH2D(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::MnvH2D*)
   {
      ::PlotUtils::MnvH2D *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PlotUtils::MnvH2D >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::MnvH2D", ::PlotUtils::MnvH2D::Class_Version(), "MnvH2D.h", 20,
                  typeid(::PlotUtils::MnvH2D), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PlotUtils::MnvH2D::Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::MnvH2D) );
      instance.SetNew(&new_PlotUtilscLcLMnvH2D);
      instance.SetNewArray(&newArray_PlotUtilscLcLMnvH2D);
      instance.SetDelete(&delete_PlotUtilscLcLMnvH2D);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLMnvH2D);
      instance.SetDestructor(&destruct_PlotUtilscLcLMnvH2D);
      instance.SetDirectoryAutoAdd(&directoryAutoAdd_PlotUtilscLcLMnvH2D);
      instance.SetMerge(&merge_PlotUtilscLcLMnvH2D);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::MnvH2D*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::MnvH2D*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::MnvH2D*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PlotUtilscLcLMnvVertErrorBand3D(void *p = nullptr);
   static void *newArray_PlotUtilscLcLMnvVertErrorBand3D(Long_t size, void *p);
   static void delete_PlotUtilscLcLMnvVertErrorBand3D(void *p);
   static void deleteArray_PlotUtilscLcLMnvVertErrorBand3D(void *p);
   static void destruct_PlotUtilscLcLMnvVertErrorBand3D(void *p);
   static void directoryAutoAdd_PlotUtilscLcLMnvVertErrorBand3D(void *obj, TDirectory *dir);
   static Long64_t merge_PlotUtilscLcLMnvVertErrorBand3D(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::MnvVertErrorBand3D*)
   {
      ::PlotUtils::MnvVertErrorBand3D *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PlotUtils::MnvVertErrorBand3D >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::MnvVertErrorBand3D", ::PlotUtils::MnvVertErrorBand3D::Class_Version(), "MnvVertErrorBand3D.h", 21,
                  typeid(::PlotUtils::MnvVertErrorBand3D), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PlotUtils::MnvVertErrorBand3D::Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::MnvVertErrorBand3D) );
      instance.SetNew(&new_PlotUtilscLcLMnvVertErrorBand3D);
      instance.SetNewArray(&newArray_PlotUtilscLcLMnvVertErrorBand3D);
      instance.SetDelete(&delete_PlotUtilscLcLMnvVertErrorBand3D);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLMnvVertErrorBand3D);
      instance.SetDestructor(&destruct_PlotUtilscLcLMnvVertErrorBand3D);
      instance.SetDirectoryAutoAdd(&directoryAutoAdd_PlotUtilscLcLMnvVertErrorBand3D);
      instance.SetMerge(&merge_PlotUtilscLcLMnvVertErrorBand3D);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::MnvVertErrorBand3D*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::MnvVertErrorBand3D*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::MnvVertErrorBand3D*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PlotUtilscLcLMnvLatErrorBand3D(void *p = nullptr);
   static void *newArray_PlotUtilscLcLMnvLatErrorBand3D(Long_t size, void *p);
   static void delete_PlotUtilscLcLMnvLatErrorBand3D(void *p);
   static void deleteArray_PlotUtilscLcLMnvLatErrorBand3D(void *p);
   static void destruct_PlotUtilscLcLMnvLatErrorBand3D(void *p);
   static void directoryAutoAdd_PlotUtilscLcLMnvLatErrorBand3D(void *obj, TDirectory *dir);
   static Long64_t merge_PlotUtilscLcLMnvLatErrorBand3D(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::MnvLatErrorBand3D*)
   {
      ::PlotUtils::MnvLatErrorBand3D *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PlotUtils::MnvLatErrorBand3D >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::MnvLatErrorBand3D", ::PlotUtils::MnvLatErrorBand3D::Class_Version(), "MnvLatErrorBand3D.h", 21,
                  typeid(::PlotUtils::MnvLatErrorBand3D), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PlotUtils::MnvLatErrorBand3D::Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::MnvLatErrorBand3D) );
      instance.SetNew(&new_PlotUtilscLcLMnvLatErrorBand3D);
      instance.SetNewArray(&newArray_PlotUtilscLcLMnvLatErrorBand3D);
      instance.SetDelete(&delete_PlotUtilscLcLMnvLatErrorBand3D);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLMnvLatErrorBand3D);
      instance.SetDestructor(&destruct_PlotUtilscLcLMnvLatErrorBand3D);
      instance.SetDirectoryAutoAdd(&directoryAutoAdd_PlotUtilscLcLMnvLatErrorBand3D);
      instance.SetMerge(&merge_PlotUtilscLcLMnvLatErrorBand3D);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::MnvLatErrorBand3D*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::MnvLatErrorBand3D*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::MnvLatErrorBand3D*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PlotUtilscLcLMnvH3D(void *p = nullptr);
   static void *newArray_PlotUtilscLcLMnvH3D(Long_t size, void *p);
   static void delete_PlotUtilscLcLMnvH3D(void *p);
   static void deleteArray_PlotUtilscLcLMnvH3D(void *p);
   static void destruct_PlotUtilscLcLMnvH3D(void *p);
   static void directoryAutoAdd_PlotUtilscLcLMnvH3D(void *obj, TDirectory *dir);
   static Long64_t merge_PlotUtilscLcLMnvH3D(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::MnvH3D*)
   {
      ::PlotUtils::MnvH3D *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PlotUtils::MnvH3D >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::MnvH3D", ::PlotUtils::MnvH3D::Class_Version(), "MnvH3D.h", 20,
                  typeid(::PlotUtils::MnvH3D), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PlotUtils::MnvH3D::Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::MnvH3D) );
      instance.SetNew(&new_PlotUtilscLcLMnvH3D);
      instance.SetNewArray(&newArray_PlotUtilscLcLMnvH3D);
      instance.SetDelete(&delete_PlotUtilscLcLMnvH3D);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLMnvH3D);
      instance.SetDestructor(&destruct_PlotUtilscLcLMnvH3D);
      instance.SetDirectoryAutoAdd(&directoryAutoAdd_PlotUtilscLcLMnvH3D);
      instance.SetMerge(&merge_PlotUtilscLcLMnvH3D);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::MnvH3D*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::MnvH3D*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::MnvH3D*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *PlotUtilscLcLMnvPlotter_Dictionary();
   static void PlotUtilscLcLMnvPlotter_TClassManip(TClass*);
   static void *new_PlotUtilscLcLMnvPlotter(void *p = nullptr);
   static void *newArray_PlotUtilscLcLMnvPlotter(Long_t size, void *p);
   static void delete_PlotUtilscLcLMnvPlotter(void *p);
   static void deleteArray_PlotUtilscLcLMnvPlotter(void *p);
   static void destruct_PlotUtilscLcLMnvPlotter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::MnvPlotter*)
   {
      ::PlotUtils::MnvPlotter *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::PlotUtils::MnvPlotter));
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::MnvPlotter", "MnvPlotter.h", 72,
                  typeid(::PlotUtils::MnvPlotter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &PlotUtilscLcLMnvPlotter_Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::MnvPlotter) );
      instance.SetNew(&new_PlotUtilscLcLMnvPlotter);
      instance.SetNewArray(&newArray_PlotUtilscLcLMnvPlotter);
      instance.SetDelete(&delete_PlotUtilscLcLMnvPlotter);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLMnvPlotter);
      instance.SetDestructor(&destruct_PlotUtilscLcLMnvPlotter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::MnvPlotter*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::MnvPlotter*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::MnvPlotter*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *PlotUtilscLcLMnvPlotter_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::PlotUtils::MnvPlotter*>(nullptr))->GetClass();
      PlotUtilscLcLMnvPlotter_TClassManip(theClass);
   return theClass;
   }

   static void PlotUtilscLcLMnvPlotter_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_PlotUtilscLcLTreeWrapper(void *p = nullptr);
   static void *newArray_PlotUtilscLcLTreeWrapper(Long_t size, void *p);
   static void delete_PlotUtilscLcLTreeWrapper(void *p);
   static void deleteArray_PlotUtilscLcLTreeWrapper(void *p);
   static void destruct_PlotUtilscLcLTreeWrapper(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::TreeWrapper*)
   {
      ::PlotUtils::TreeWrapper *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PlotUtils::TreeWrapper >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::TreeWrapper", ::PlotUtils::TreeWrapper::Class_Version(), "TreeWrapper.h", 67,
                  typeid(::PlotUtils::TreeWrapper), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PlotUtils::TreeWrapper::Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::TreeWrapper) );
      instance.SetNew(&new_PlotUtilscLcLTreeWrapper);
      instance.SetNewArray(&newArray_PlotUtilscLcLTreeWrapper);
      instance.SetDelete(&delete_PlotUtilscLcLTreeWrapper);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLTreeWrapper);
      instance.SetDestructor(&destruct_PlotUtilscLcLTreeWrapper);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::TreeWrapper*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::TreeWrapper*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::TreeWrapper*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_PlotUtilscLcLChainWrapper(void *p);
   static void deleteArray_PlotUtilscLcLChainWrapper(void *p);
   static void destruct_PlotUtilscLcLChainWrapper(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::ChainWrapper*)
   {
      ::PlotUtils::ChainWrapper *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PlotUtils::ChainWrapper >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::ChainWrapper", ::PlotUtils::ChainWrapper::Class_Version(), "ChainWrapper.h", 13,
                  typeid(::PlotUtils::ChainWrapper), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PlotUtils::ChainWrapper::Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::ChainWrapper) );
      instance.SetDelete(&delete_PlotUtilscLcLChainWrapper);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLChainWrapper);
      instance.SetDestructor(&destruct_PlotUtilscLcLChainWrapper);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::ChainWrapper*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::ChainWrapper*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::ChainWrapper*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PlotUtilscLcLGridCanvas(void *p = nullptr);
   static void *newArray_PlotUtilscLcLGridCanvas(Long_t size, void *p);
   static void delete_PlotUtilscLcLGridCanvas(void *p);
   static void deleteArray_PlotUtilscLcLGridCanvas(void *p);
   static void destruct_PlotUtilscLcLGridCanvas(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::GridCanvas*)
   {
      ::PlotUtils::GridCanvas *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PlotUtils::GridCanvas >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::GridCanvas", ::PlotUtils::GridCanvas::Class_Version(), "GridCanvas.h", 17,
                  typeid(::PlotUtils::GridCanvas), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PlotUtils::GridCanvas::Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::GridCanvas) );
      instance.SetNew(&new_PlotUtilscLcLGridCanvas);
      instance.SetNewArray(&newArray_PlotUtilscLcLGridCanvas);
      instance.SetDelete(&delete_PlotUtilscLcLGridCanvas);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLGridCanvas);
      instance.SetDestructor(&destruct_PlotUtilscLcLGridCanvas);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::GridCanvas*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::GridCanvas*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::GridCanvas*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *PlotUtilscLcLFluxReweighter_Dictionary();
   static void PlotUtilscLcLFluxReweighter_TClassManip(TClass*);
   static void delete_PlotUtilscLcLFluxReweighter(void *p);
   static void deleteArray_PlotUtilscLcLFluxReweighter(void *p);
   static void destruct_PlotUtilscLcLFluxReweighter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::FluxReweighter*)
   {
      ::PlotUtils::FluxReweighter *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::PlotUtils::FluxReweighter));
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::FluxReweighter", "FluxReweighter.h", 68,
                  typeid(::PlotUtils::FluxReweighter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &PlotUtilscLcLFluxReweighter_Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::FluxReweighter) );
      instance.SetDelete(&delete_PlotUtilscLcLFluxReweighter);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLFluxReweighter);
      instance.SetDestructor(&destruct_PlotUtilscLcLFluxReweighter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::FluxReweighter*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::FluxReweighter*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::FluxReweighter*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *PlotUtilscLcLFluxReweighter_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::PlotUtils::FluxReweighter*>(nullptr))->GetClass();
      PlotUtilscLcLFluxReweighter_TClassManip(theClass);
   return theClass;
   }

   static void PlotUtilscLcLFluxReweighter_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *PlotUtilscLcLBaseUniverse_Dictionary();
   static void PlotUtilscLcLBaseUniverse_TClassManip(TClass*);
   static void *new_PlotUtilscLcLBaseUniverse(void *p = nullptr);
   static void *newArray_PlotUtilscLcLBaseUniverse(Long_t size, void *p);
   static void delete_PlotUtilscLcLBaseUniverse(void *p);
   static void deleteArray_PlotUtilscLcLBaseUniverse(void *p);
   static void destruct_PlotUtilscLcLBaseUniverse(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotUtils::BaseUniverse*)
   {
      ::PlotUtils::BaseUniverse *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::PlotUtils::BaseUniverse));
      static ::ROOT::TGenericClassInfo 
         instance("PlotUtils::BaseUniverse", "BaseUniverse.h", 10,
                  typeid(::PlotUtils::BaseUniverse), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &PlotUtilscLcLBaseUniverse_Dictionary, isa_proxy, 4,
                  sizeof(::PlotUtils::BaseUniverse) );
      instance.SetNew(&new_PlotUtilscLcLBaseUniverse);
      instance.SetNewArray(&newArray_PlotUtilscLcLBaseUniverse);
      instance.SetDelete(&delete_PlotUtilscLcLBaseUniverse);
      instance.SetDeleteArray(&deleteArray_PlotUtilscLcLBaseUniverse);
      instance.SetDestructor(&destruct_PlotUtilscLcLBaseUniverse);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotUtils::BaseUniverse*)
   {
      return GenerateInitInstanceLocal(static_cast<::PlotUtils::BaseUniverse*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PlotUtils::BaseUniverse*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *PlotUtilscLcLBaseUniverse_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::PlotUtils::BaseUniverse*>(nullptr))->GetClass();
      PlotUtilscLcLBaseUniverse_TClassManip(theClass);
   return theClass;
   }

   static void PlotUtilscLcLBaseUniverse_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ROOTcLcLexception_Dictionary();
   static void ROOTcLcLexception_TClassManip(TClass*);
   static void delete_ROOTcLcLexception(void *p);
   static void deleteArray_ROOTcLcLexception(void *p);
   static void destruct_ROOTcLcLexception(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ROOT::exception*)
   {
      ::ROOT::exception *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ROOT::exception));
      static ::ROOT::TGenericClassInfo 
         instance("ROOT::exception", "ErrorHandler.h", 50,
                  typeid(::ROOT::exception), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ROOTcLcLexception_Dictionary, isa_proxy, 4,
                  sizeof(::ROOT::exception) );
      instance.SetDelete(&delete_ROOTcLcLexception);
      instance.SetDeleteArray(&deleteArray_ROOTcLcLexception);
      instance.SetDestructor(&destruct_ROOTcLcLexception);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ROOT::exception*)
   {
      return GenerateInitInstanceLocal(static_cast<::ROOT::exception*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ROOT::exception*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ROOTcLcLexception_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::ROOT::exception*>(nullptr))->GetClass();
      ROOTcLcLexception_TClassManip(theClass);
   return theClass;
   }

   static void ROOTcLcLexception_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ROOTcLcLerror_Dictionary();
   static void ROOTcLcLerror_TClassManip(TClass*);
   static void delete_ROOTcLcLerror(void *p);
   static void deleteArray_ROOTcLcLerror(void *p);
   static void destruct_ROOTcLcLerror(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ROOT::error*)
   {
      ::ROOT::error *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ROOT::error));
      static ::ROOT::TGenericClassInfo 
         instance("ROOT::error", "ErrorHandler.h", 63,
                  typeid(::ROOT::error), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ROOTcLcLerror_Dictionary, isa_proxy, 4,
                  sizeof(::ROOT::error) );
      instance.SetDelete(&delete_ROOTcLcLerror);
      instance.SetDeleteArray(&deleteArray_ROOTcLcLerror);
      instance.SetDestructor(&destruct_ROOTcLcLerror);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ROOT::error*)
   {
      return GenerateInitInstanceLocal(static_cast<::ROOT::error*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ROOT::error*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ROOTcLcLerror_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::ROOT::error*>(nullptr))->GetClass();
      ROOTcLcLerror_TClassManip(theClass);
   return theClass;
   }

   static void ROOTcLcLerror_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ROOTcLcLwarning_Dictionary();
   static void ROOTcLcLwarning_TClassManip(TClass*);
   static void delete_ROOTcLcLwarning(void *p);
   static void deleteArray_ROOTcLcLwarning(void *p);
   static void destruct_ROOTcLcLwarning(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ROOT::warning*)
   {
      ::ROOT::warning *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ROOT::warning));
      static ::ROOT::TGenericClassInfo 
         instance("ROOT::warning", "ErrorHandler.h", 72,
                  typeid(::ROOT::warning), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ROOTcLcLwarning_Dictionary, isa_proxy, 4,
                  sizeof(::ROOT::warning) );
      instance.SetDelete(&delete_ROOTcLcLwarning);
      instance.SetDeleteArray(&deleteArray_ROOTcLcLwarning);
      instance.SetDestructor(&destruct_ROOTcLcLwarning);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ROOT::warning*)
   {
      return GenerateInitInstanceLocal(static_cast<::ROOT::warning*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ROOT::warning*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ROOTcLcLwarning_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::ROOT::warning*>(nullptr))->GetClass();
      ROOTcLcLwarning_TClassManip(theClass);
   return theClass;
   }

   static void ROOTcLcLwarning_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ROOTcLcLdetailcLcLbeforeMain_Dictionary();
   static void ROOTcLcLdetailcLcLbeforeMain_TClassManip(TClass*);
   static void *new_ROOTcLcLdetailcLcLbeforeMain(void *p = nullptr);
   static void *newArray_ROOTcLcLdetailcLcLbeforeMain(Long_t size, void *p);
   static void delete_ROOTcLcLdetailcLcLbeforeMain(void *p);
   static void deleteArray_ROOTcLcLdetailcLcLbeforeMain(void *p);
   static void destruct_ROOTcLcLdetailcLcLbeforeMain(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ROOT::detail::beforeMain*)
   {
      ::ROOT::detail::beforeMain *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ROOT::detail::beforeMain));
      static ::ROOT::TGenericClassInfo 
         instance("ROOT::detail::beforeMain", "ErrorHandler.h", 84,
                  typeid(::ROOT::detail::beforeMain), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ROOTcLcLdetailcLcLbeforeMain_Dictionary, isa_proxy, 4,
                  sizeof(::ROOT::detail::beforeMain) );
      instance.SetNew(&new_ROOTcLcLdetailcLcLbeforeMain);
      instance.SetNewArray(&newArray_ROOTcLcLdetailcLcLbeforeMain);
      instance.SetDelete(&delete_ROOTcLcLdetailcLcLbeforeMain);
      instance.SetDeleteArray(&deleteArray_ROOTcLcLdetailcLcLbeforeMain);
      instance.SetDestructor(&destruct_ROOTcLcLdetailcLcLbeforeMain);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ROOT::detail::beforeMain*)
   {
      return GenerateInitInstanceLocal(static_cast<::ROOT::detail::beforeMain*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ROOT::detail::beforeMain*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ROOTcLcLdetailcLcLbeforeMain_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::ROOT::detail::beforeMain*>(nullptr))->GetClass();
      ROOTcLcLdetailcLcLbeforeMain_TClassManip(theClass);
   return theClass;
   }

   static void ROOTcLcLdetailcLcLbeforeMain_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *pairlEstringcOvectorlEstringgRsPgR_Dictionary();
   static void pairlEstringcOvectorlEstringgRsPgR_TClassManip(TClass*);
   static void *new_pairlEstringcOvectorlEstringgRsPgR(void *p = nullptr);
   static void *newArray_pairlEstringcOvectorlEstringgRsPgR(Long_t size, void *p);
   static void delete_pairlEstringcOvectorlEstringgRsPgR(void *p);
   static void deleteArray_pairlEstringcOvectorlEstringgRsPgR(void *p);
   static void destruct_pairlEstringcOvectorlEstringgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const pair<string,vector<string> >*)
   {
      pair<string,vector<string> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(pair<string,vector<string> >));
      static ::ROOT::TGenericClassInfo 
         instance("pair<string,vector<string> >", "", 72,
                  typeid(pair<string,vector<string> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &pairlEstringcOvectorlEstringgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(pair<string,vector<string> >) );
      instance.SetNew(&new_pairlEstringcOvectorlEstringgRsPgR);
      instance.SetNewArray(&newArray_pairlEstringcOvectorlEstringgRsPgR);
      instance.SetDelete(&delete_pairlEstringcOvectorlEstringgRsPgR);
      instance.SetDeleteArray(&deleteArray_pairlEstringcOvectorlEstringgRsPgR);
      instance.SetDestructor(&destruct_pairlEstringcOvectorlEstringgRsPgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("pair<string,vector<string> >","pair<std::string,std::vector<std::string> >"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("pair<string,vector<string> >","std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const pair<string,vector<string> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *pairlEstringcOvectorlEstringgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const pair<string,vector<string> >*>(nullptr))->GetClass();
      pairlEstringcOvectorlEstringgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void pairlEstringcOvectorlEstringgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR_Dictionary();
   static void pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR_TClassManip(TClass*);
   static void *new_pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(void *p = nullptr);
   static void *newArray_pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(Long_t size, void *p);
   static void delete_pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(void *p);
   static void deleteArray_pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(void *p);
   static void destruct_pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const pair<string,PlotUtils::MnvLatErrorBand*>*)
   {
      pair<string,PlotUtils::MnvLatErrorBand*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(pair<string,PlotUtils::MnvLatErrorBand*>));
      static ::ROOT::TGenericClassInfo 
         instance("pair<string,PlotUtils::MnvLatErrorBand*>", "", 75,
                  typeid(pair<string,PlotUtils::MnvLatErrorBand*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR_Dictionary, isa_proxy, 4,
                  sizeof(pair<string,PlotUtils::MnvLatErrorBand*>) );
      instance.SetNew(&new_pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR);
      instance.SetNewArray(&newArray_pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR);
      instance.SetDelete(&delete_pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR);
      instance.SetDeleteArray(&deleteArray_pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR);
      instance.SetDestructor(&destruct_pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("pair<string,PlotUtils::MnvLatErrorBand*>","pair<std::string,PlotUtils::MnvLatErrorBand*>"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("pair<string,PlotUtils::MnvLatErrorBand*>","std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, PlotUtils::MnvLatErrorBand*>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const pair<string,PlotUtils::MnvLatErrorBand*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const pair<string,PlotUtils::MnvLatErrorBand*>*>(nullptr))->GetClass();
      pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR_TClassManip(theClass);
   return theClass;
   }

   static void pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR_Dictionary();
   static void pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR_TClassManip(TClass*);
   static void *new_pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(void *p = nullptr);
   static void *newArray_pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(Long_t size, void *p);
   static void delete_pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(void *p);
   static void deleteArray_pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(void *p);
   static void destruct_pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const pair<string,PlotUtils::MnvVertErrorBand*>*)
   {
      pair<string,PlotUtils::MnvVertErrorBand*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(pair<string,PlotUtils::MnvVertErrorBand*>));
      static ::ROOT::TGenericClassInfo 
         instance("pair<string,PlotUtils::MnvVertErrorBand*>", "", 76,
                  typeid(pair<string,PlotUtils::MnvVertErrorBand*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR_Dictionary, isa_proxy, 4,
                  sizeof(pair<string,PlotUtils::MnvVertErrorBand*>) );
      instance.SetNew(&new_pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR);
      instance.SetNewArray(&newArray_pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR);
      instance.SetDelete(&delete_pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR);
      instance.SetDeleteArray(&deleteArray_pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR);
      instance.SetDestructor(&destruct_pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("pair<string,PlotUtils::MnvVertErrorBand*>","pair<std::string,PlotUtils::MnvVertErrorBand*>"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("pair<string,PlotUtils::MnvVertErrorBand*>","std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, PlotUtils::MnvVertErrorBand*>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const pair<string,PlotUtils::MnvVertErrorBand*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const pair<string,PlotUtils::MnvVertErrorBand*>*>(nullptr))->GetClass();
      pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR_TClassManip(theClass);
   return theClass;
   }

   static void pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR_Dictionary();
   static void pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR_TClassManip(TClass*);
   static void *new_pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(void *p = nullptr);
   static void *newArray_pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(Long_t size, void *p);
   static void delete_pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(void *p);
   static void deleteArray_pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(void *p);
   static void destruct_pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const pair<string,PlotUtils::MnvVertErrorBand3D*>*)
   {
      pair<string,PlotUtils::MnvVertErrorBand3D*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(pair<string,PlotUtils::MnvVertErrorBand3D*>));
      static ::ROOT::TGenericClassInfo 
         instance("pair<string,PlotUtils::MnvVertErrorBand3D*>", "", 77,
                  typeid(pair<string,PlotUtils::MnvVertErrorBand3D*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR_Dictionary, isa_proxy, 4,
                  sizeof(pair<string,PlotUtils::MnvVertErrorBand3D*>) );
      instance.SetNew(&new_pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR);
      instance.SetNewArray(&newArray_pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR);
      instance.SetDelete(&delete_pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR);
      instance.SetDeleteArray(&deleteArray_pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR);
      instance.SetDestructor(&destruct_pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("pair<string,PlotUtils::MnvVertErrorBand3D*>","pair<std::string,PlotUtils::MnvVertErrorBand3D*>"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("pair<string,PlotUtils::MnvVertErrorBand3D*>","std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, PlotUtils::MnvVertErrorBand3D*>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const pair<string,PlotUtils::MnvVertErrorBand3D*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const pair<string,PlotUtils::MnvVertErrorBand3D*>*>(nullptr))->GetClass();
      pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR_TClassManip(theClass);
   return theClass;
   }

   static void pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR_Dictionary();
   static void pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR_TClassManip(TClass*);
   static void *new_pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(void *p = nullptr);
   static void *newArray_pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(Long_t size, void *p);
   static void delete_pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(void *p);
   static void deleteArray_pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(void *p);
   static void destruct_pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const pair<string,PlotUtils::MnvLatErrorBand3D*>*)
   {
      pair<string,PlotUtils::MnvLatErrorBand3D*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(pair<string,PlotUtils::MnvLatErrorBand3D*>));
      static ::ROOT::TGenericClassInfo 
         instance("pair<string,PlotUtils::MnvLatErrorBand3D*>", "", 78,
                  typeid(pair<string,PlotUtils::MnvLatErrorBand3D*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR_Dictionary, isa_proxy, 4,
                  sizeof(pair<string,PlotUtils::MnvLatErrorBand3D*>) );
      instance.SetNew(&new_pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR);
      instance.SetNewArray(&newArray_pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR);
      instance.SetDelete(&delete_pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR);
      instance.SetDeleteArray(&deleteArray_pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR);
      instance.SetDestructor(&destruct_pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("pair<string,PlotUtils::MnvLatErrorBand3D*>","pair<std::string,PlotUtils::MnvLatErrorBand3D*>"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("pair<string,PlotUtils::MnvLatErrorBand3D*>","std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, PlotUtils::MnvLatErrorBand3D*>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const pair<string,PlotUtils::MnvLatErrorBand3D*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const pair<string,PlotUtils::MnvLatErrorBand3D*>*>(nullptr))->GetClass();
      pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR_TClassManip(theClass);
   return theClass;
   }

   static void pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *pairlEstringcOTMatrixTlEdoublegRmUgR_Dictionary();
   static void pairlEstringcOTMatrixTlEdoublegRmUgR_TClassManip(TClass*);
   static void *new_pairlEstringcOTMatrixTlEdoublegRmUgR(void *p = nullptr);
   static void *newArray_pairlEstringcOTMatrixTlEdoublegRmUgR(Long_t size, void *p);
   static void delete_pairlEstringcOTMatrixTlEdoublegRmUgR(void *p);
   static void deleteArray_pairlEstringcOTMatrixTlEdoublegRmUgR(void *p);
   static void destruct_pairlEstringcOTMatrixTlEdoublegRmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const pair<string,TMatrixT<double>*>*)
   {
      pair<string,TMatrixT<double>*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(pair<string,TMatrixT<double>*>));
      static ::ROOT::TGenericClassInfo 
         instance("pair<string,TMatrixT<double>*>", "", 79,
                  typeid(pair<string,TMatrixT<double>*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &pairlEstringcOTMatrixTlEdoublegRmUgR_Dictionary, isa_proxy, 4,
                  sizeof(pair<string,TMatrixT<double>*>) );
      instance.SetNew(&new_pairlEstringcOTMatrixTlEdoublegRmUgR);
      instance.SetNewArray(&newArray_pairlEstringcOTMatrixTlEdoublegRmUgR);
      instance.SetDelete(&delete_pairlEstringcOTMatrixTlEdoublegRmUgR);
      instance.SetDeleteArray(&deleteArray_pairlEstringcOTMatrixTlEdoublegRmUgR);
      instance.SetDestructor(&destruct_pairlEstringcOTMatrixTlEdoublegRmUgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("pair<string,TMatrixT<double>*>","pair<std::string,TMatrixT<double>*>"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("pair<string,TMatrixT<double>*>","std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, TMatrixT<double>*>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const pair<string,TMatrixT<double>*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *pairlEstringcOTMatrixTlEdoublegRmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const pair<string,TMatrixT<double>*>*>(nullptr))->GetClass();
      pairlEstringcOTMatrixTlEdoublegRmUgR_TClassManip(theClass);
   return theClass;
   }

   static void pairlEstringcOTMatrixTlEdoublegRmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR_Dictionary();
   static void pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR_TClassManip(TClass*);
   static void *new_pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(void *p = nullptr);
   static void *newArray_pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(Long_t size, void *p);
   static void delete_pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(void *p);
   static void deleteArray_pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(void *p);
   static void destruct_pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const pair<string,PlotUtils::MnvVertErrorBand2D*>*)
   {
      pair<string,PlotUtils::MnvVertErrorBand2D*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(pair<string,PlotUtils::MnvVertErrorBand2D*>));
      static ::ROOT::TGenericClassInfo 
         instance("pair<string,PlotUtils::MnvVertErrorBand2D*>", "", 80,
                  typeid(pair<string,PlotUtils::MnvVertErrorBand2D*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR_Dictionary, isa_proxy, 4,
                  sizeof(pair<string,PlotUtils::MnvVertErrorBand2D*>) );
      instance.SetNew(&new_pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR);
      instance.SetNewArray(&newArray_pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR);
      instance.SetDelete(&delete_pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR);
      instance.SetDeleteArray(&deleteArray_pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR);
      instance.SetDestructor(&destruct_pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("pair<string,PlotUtils::MnvVertErrorBand2D*>","pair<std::string,PlotUtils::MnvVertErrorBand2D*>"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("pair<string,PlotUtils::MnvVertErrorBand2D*>","std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, PlotUtils::MnvVertErrorBand2D*>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const pair<string,PlotUtils::MnvVertErrorBand2D*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const pair<string,PlotUtils::MnvVertErrorBand2D*>*>(nullptr))->GetClass();
      pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR_TClassManip(theClass);
   return theClass;
   }

   static void pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR_Dictionary();
   static void pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR_TClassManip(TClass*);
   static void *new_pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(void *p = nullptr);
   static void *newArray_pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(Long_t size, void *p);
   static void delete_pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(void *p);
   static void deleteArray_pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(void *p);
   static void destruct_pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const pair<string,PlotUtils::MnvLatErrorBand2D*>*)
   {
      pair<string,PlotUtils::MnvLatErrorBand2D*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(pair<string,PlotUtils::MnvLatErrorBand2D*>));
      static ::ROOT::TGenericClassInfo 
         instance("pair<string,PlotUtils::MnvLatErrorBand2D*>", "", 81,
                  typeid(pair<string,PlotUtils::MnvLatErrorBand2D*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR_Dictionary, isa_proxy, 4,
                  sizeof(pair<string,PlotUtils::MnvLatErrorBand2D*>) );
      instance.SetNew(&new_pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR);
      instance.SetNewArray(&newArray_pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR);
      instance.SetDelete(&delete_pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR);
      instance.SetDeleteArray(&deleteArray_pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR);
      instance.SetDestructor(&destruct_pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("pair<string,PlotUtils::MnvLatErrorBand2D*>","pair<std::string,PlotUtils::MnvLatErrorBand2D*>"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("pair<string,PlotUtils::MnvLatErrorBand2D*>","std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, PlotUtils::MnvLatErrorBand2D*>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const pair<string,PlotUtils::MnvLatErrorBand2D*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const pair<string,PlotUtils::MnvLatErrorBand2D*>*>(nullptr))->GetClass();
      pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR_TClassManip(theClass);
   return theClass;
   }

   static void pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace PlotUtils {
//______________________________________________________________________________
atomic_TClass_ptr MnvVertErrorBand::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *MnvVertErrorBand::Class_Name()
{
   return "PlotUtils::MnvVertErrorBand";
}

//______________________________________________________________________________
const char *MnvVertErrorBand::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvVertErrorBand*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int MnvVertErrorBand::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvVertErrorBand*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MnvVertErrorBand::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvVertErrorBand*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MnvVertErrorBand::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvVertErrorBand*)nullptr)->GetClass(); }
   return fgIsA;
}

} // namespace PlotUtils
namespace PlotUtils {
//______________________________________________________________________________
atomic_TClass_ptr MnvLatErrorBand::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *MnvLatErrorBand::Class_Name()
{
   return "PlotUtils::MnvLatErrorBand";
}

//______________________________________________________________________________
const char *MnvLatErrorBand::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvLatErrorBand*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int MnvLatErrorBand::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvLatErrorBand*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MnvLatErrorBand::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvLatErrorBand*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MnvLatErrorBand::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvLatErrorBand*)nullptr)->GetClass(); }
   return fgIsA;
}

} // namespace PlotUtils
namespace PlotUtils {
//______________________________________________________________________________
atomic_TClass_ptr MnvH1D::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *MnvH1D::Class_Name()
{
   return "PlotUtils::MnvH1D";
}

//______________________________________________________________________________
const char *MnvH1D::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvH1D*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int MnvH1D::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvH1D*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MnvH1D::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvH1D*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MnvH1D::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvH1D*)nullptr)->GetClass(); }
   return fgIsA;
}

} // namespace PlotUtils
namespace PlotUtils {
//______________________________________________________________________________
atomic_TClass_ptr MnvVertErrorBand2D::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *MnvVertErrorBand2D::Class_Name()
{
   return "PlotUtils::MnvVertErrorBand2D";
}

//______________________________________________________________________________
const char *MnvVertErrorBand2D::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvVertErrorBand2D*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int MnvVertErrorBand2D::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvVertErrorBand2D*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MnvVertErrorBand2D::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvVertErrorBand2D*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MnvVertErrorBand2D::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvVertErrorBand2D*)nullptr)->GetClass(); }
   return fgIsA;
}

} // namespace PlotUtils
namespace PlotUtils {
//______________________________________________________________________________
atomic_TClass_ptr MnvLatErrorBand2D::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *MnvLatErrorBand2D::Class_Name()
{
   return "PlotUtils::MnvLatErrorBand2D";
}

//______________________________________________________________________________
const char *MnvLatErrorBand2D::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvLatErrorBand2D*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int MnvLatErrorBand2D::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvLatErrorBand2D*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MnvLatErrorBand2D::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvLatErrorBand2D*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MnvLatErrorBand2D::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvLatErrorBand2D*)nullptr)->GetClass(); }
   return fgIsA;
}

} // namespace PlotUtils
namespace PlotUtils {
//______________________________________________________________________________
atomic_TClass_ptr MnvH2D::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *MnvH2D::Class_Name()
{
   return "PlotUtils::MnvH2D";
}

//______________________________________________________________________________
const char *MnvH2D::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvH2D*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int MnvH2D::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvH2D*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MnvH2D::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvH2D*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MnvH2D::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvH2D*)nullptr)->GetClass(); }
   return fgIsA;
}

} // namespace PlotUtils
namespace PlotUtils {
//______________________________________________________________________________
atomic_TClass_ptr MnvVertErrorBand3D::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *MnvVertErrorBand3D::Class_Name()
{
   return "PlotUtils::MnvVertErrorBand3D";
}

//______________________________________________________________________________
const char *MnvVertErrorBand3D::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvVertErrorBand3D*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int MnvVertErrorBand3D::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvVertErrorBand3D*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MnvVertErrorBand3D::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvVertErrorBand3D*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MnvVertErrorBand3D::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvVertErrorBand3D*)nullptr)->GetClass(); }
   return fgIsA;
}

} // namespace PlotUtils
namespace PlotUtils {
//______________________________________________________________________________
atomic_TClass_ptr MnvLatErrorBand3D::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *MnvLatErrorBand3D::Class_Name()
{
   return "PlotUtils::MnvLatErrorBand3D";
}

//______________________________________________________________________________
const char *MnvLatErrorBand3D::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvLatErrorBand3D*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int MnvLatErrorBand3D::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvLatErrorBand3D*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MnvLatErrorBand3D::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvLatErrorBand3D*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MnvLatErrorBand3D::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvLatErrorBand3D*)nullptr)->GetClass(); }
   return fgIsA;
}

} // namespace PlotUtils
namespace PlotUtils {
//______________________________________________________________________________
atomic_TClass_ptr MnvH3D::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *MnvH3D::Class_Name()
{
   return "PlotUtils::MnvH3D";
}

//______________________________________________________________________________
const char *MnvH3D::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvH3D*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int MnvH3D::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvH3D*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MnvH3D::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvH3D*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MnvH3D::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::MnvH3D*)nullptr)->GetClass(); }
   return fgIsA;
}

} // namespace PlotUtils
namespace PlotUtils {
//______________________________________________________________________________
atomic_TClass_ptr TreeWrapper::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TreeWrapper::Class_Name()
{
   return "PlotUtils::TreeWrapper";
}

//______________________________________________________________________________
const char *TreeWrapper::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::TreeWrapper*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TreeWrapper::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::TreeWrapper*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TreeWrapper::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::TreeWrapper*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TreeWrapper::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::TreeWrapper*)nullptr)->GetClass(); }
   return fgIsA;
}

} // namespace PlotUtils
namespace PlotUtils {
//______________________________________________________________________________
atomic_TClass_ptr ChainWrapper::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *ChainWrapper::Class_Name()
{
   return "PlotUtils::ChainWrapper";
}

//______________________________________________________________________________
const char *ChainWrapper::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::ChainWrapper*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int ChainWrapper::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::ChainWrapper*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ChainWrapper::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::ChainWrapper*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ChainWrapper::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::ChainWrapper*)nullptr)->GetClass(); }
   return fgIsA;
}

} // namespace PlotUtils
namespace PlotUtils {
//______________________________________________________________________________
atomic_TClass_ptr GridCanvas::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *GridCanvas::Class_Name()
{
   return "PlotUtils::GridCanvas";
}

//______________________________________________________________________________
const char *GridCanvas::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::GridCanvas*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int GridCanvas::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::GridCanvas*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GridCanvas::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::GridCanvas*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GridCanvas::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotUtils::GridCanvas*)nullptr)->GetClass(); }
   return fgIsA;
}

} // namespace PlotUtils
namespace ROOT {
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLBadUnivWgtError(void *p) {
      delete (static_cast<::PlotUtils::BadUnivWgtError*>(p));
   }
   static void deleteArray_PlotUtilscLcLBadUnivWgtError(void *p) {
      delete [] (static_cast<::PlotUtils::BadUnivWgtError*>(p));
   }
   static void destruct_PlotUtilscLcLBadUnivWgtError(void *p) {
      typedef ::PlotUtils::BadUnivWgtError current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::PlotUtils::BadUnivWgtError

namespace ROOT {
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLConstraintLoadError(void *p) {
      delete (static_cast<::PlotUtils::ConstraintLoadError*>(p));
   }
   static void deleteArray_PlotUtilscLcLConstraintLoadError(void *p) {
      delete [] (static_cast<::PlotUtils::ConstraintLoadError*>(p));
   }
   static void destruct_PlotUtilscLcLConstraintLoadError(void *p) {
      typedef ::PlotUtils::ConstraintLoadError current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::PlotUtils::ConstraintLoadError

namespace ROOT {
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLConstraintAccessError(void *p) {
      delete (static_cast<::PlotUtils::ConstraintAccessError*>(p));
   }
   static void deleteArray_PlotUtilscLcLConstraintAccessError(void *p) {
      delete [] (static_cast<::PlotUtils::ConstraintAccessError*>(p));
   }
   static void destruct_PlotUtilscLcLConstraintAccessError(void *p) {
      typedef ::PlotUtils::ConstraintAccessError current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::PlotUtils::ConstraintAccessError

namespace ROOT {
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLMissingFluxUnivError(void *p) {
      delete (static_cast<::PlotUtils::MissingFluxUnivError*>(p));
   }
   static void deleteArray_PlotUtilscLcLMissingFluxUnivError(void *p) {
      delete [] (static_cast<::PlotUtils::MissingFluxUnivError*>(p));
   }
   static void destruct_PlotUtilscLcLMissingFluxUnivError(void *p) {
      typedef ::PlotUtils::MissingFluxUnivError current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::PlotUtils::MissingFluxUnivError

namespace ROOT {
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLMissingSpectatorStrategyError(void *p) {
      delete (static_cast<::PlotUtils::MissingSpectatorStrategyError*>(p));
   }
   static void deleteArray_PlotUtilscLcLMissingSpectatorStrategyError(void *p) {
      delete [] (static_cast<::PlotUtils::MissingSpectatorStrategyError*>(p));
   }
   static void destruct_PlotUtilscLcLMissingSpectatorStrategyError(void *p) {
      typedef ::PlotUtils::MissingSpectatorStrategyError current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::PlotUtils::MissingSpectatorStrategyError

namespace ROOT {
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLNoFluxUnivError(void *p) {
      delete (static_cast<::PlotUtils::NoFluxUnivError*>(p));
   }
   static void deleteArray_PlotUtilscLcLNoFluxUnivError(void *p) {
      delete [] (static_cast<::PlotUtils::NoFluxUnivError*>(p));
   }
   static void destruct_PlotUtilscLcLNoFluxUnivError(void *p) {
      typedef ::PlotUtils::NoFluxUnivError current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::PlotUtils::NoFluxUnivError

namespace ROOT {
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLNoWgtdSpreadError(void *p) {
      delete (static_cast<::PlotUtils::NoWgtdSpreadError*>(p));
   }
   static void deleteArray_PlotUtilscLcLNoWgtdSpreadError(void *p) {
      delete [] (static_cast<::PlotUtils::NoWgtdSpreadError*>(p));
   }
   static void destruct_PlotUtilscLcLNoWgtdSpreadError(void *p) {
      typedef ::PlotUtils::NoWgtdSpreadError current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::PlotUtils::NoWgtdSpreadError

namespace ROOT {
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLSpectatorConstraintCollisionError(void *p) {
      delete (static_cast<::PlotUtils::SpectatorConstraintCollisionError*>(p));
   }
   static void deleteArray_PlotUtilscLcLSpectatorConstraintCollisionError(void *p) {
      delete [] (static_cast<::PlotUtils::SpectatorConstraintCollisionError*>(p));
   }
   static void destruct_PlotUtilscLcLSpectatorConstraintCollisionError(void *p) {
      typedef ::PlotUtils::SpectatorConstraintCollisionError current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::PlotUtils::SpectatorConstraintCollisionError

namespace PlotUtils {
//______________________________________________________________________________
void MnvVertErrorBand::Streamer(TBuffer &R__b)
{
   // Stream an object of class PlotUtils::MnvVertErrorBand.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PlotUtils::MnvVertErrorBand::Class(),this);
   } else {
      R__b.WriteClassBuffer(PlotUtils::MnvVertErrorBand::Class(),this);
   }
}

} // namespace PlotUtils
namespace ROOT {
   // Wrappers around operator new
   static void *new_PlotUtilscLcLMnvVertErrorBand(void *p) {
      return  p ? new(p) ::PlotUtils::MnvVertErrorBand : new ::PlotUtils::MnvVertErrorBand;
   }
   static void *newArray_PlotUtilscLcLMnvVertErrorBand(Long_t nElements, void *p) {
      return p ? new(p) ::PlotUtils::MnvVertErrorBand[nElements] : new ::PlotUtils::MnvVertErrorBand[nElements];
   }
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLMnvVertErrorBand(void *p) {
      delete (static_cast<::PlotUtils::MnvVertErrorBand*>(p));
   }
   static void deleteArray_PlotUtilscLcLMnvVertErrorBand(void *p) {
      delete [] (static_cast<::PlotUtils::MnvVertErrorBand*>(p));
   }
   static void destruct_PlotUtilscLcLMnvVertErrorBand(void *p) {
      typedef ::PlotUtils::MnvVertErrorBand current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
   // Wrapper around the directory auto add.
   static void directoryAutoAdd_PlotUtilscLcLMnvVertErrorBand(void *p, TDirectory *dir) {
      ((::PlotUtils::MnvVertErrorBand*)p)->DirectoryAutoAdd(dir);
   }
   // Wrapper around the merge function.
   static Long64_t  merge_PlotUtilscLcLMnvVertErrorBand(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::PlotUtils::MnvVertErrorBand*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::PlotUtils::MnvVertErrorBand

namespace PlotUtils {
//______________________________________________________________________________
void MnvLatErrorBand::Streamer(TBuffer &R__b)
{
   // Stream an object of class PlotUtils::MnvLatErrorBand.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PlotUtils::MnvLatErrorBand::Class(),this);
   } else {
      R__b.WriteClassBuffer(PlotUtils::MnvLatErrorBand::Class(),this);
   }
}

} // namespace PlotUtils
namespace ROOT {
   // Wrappers around operator new
   static void *new_PlotUtilscLcLMnvLatErrorBand(void *p) {
      return  p ? new(p) ::PlotUtils::MnvLatErrorBand : new ::PlotUtils::MnvLatErrorBand;
   }
   static void *newArray_PlotUtilscLcLMnvLatErrorBand(Long_t nElements, void *p) {
      return p ? new(p) ::PlotUtils::MnvLatErrorBand[nElements] : new ::PlotUtils::MnvLatErrorBand[nElements];
   }
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLMnvLatErrorBand(void *p) {
      delete (static_cast<::PlotUtils::MnvLatErrorBand*>(p));
   }
   static void deleteArray_PlotUtilscLcLMnvLatErrorBand(void *p) {
      delete [] (static_cast<::PlotUtils::MnvLatErrorBand*>(p));
   }
   static void destruct_PlotUtilscLcLMnvLatErrorBand(void *p) {
      typedef ::PlotUtils::MnvLatErrorBand current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
   // Wrapper around the directory auto add.
   static void directoryAutoAdd_PlotUtilscLcLMnvLatErrorBand(void *p, TDirectory *dir) {
      ((::PlotUtils::MnvLatErrorBand*)p)->DirectoryAutoAdd(dir);
   }
   // Wrapper around the merge function.
   static Long64_t  merge_PlotUtilscLcLMnvLatErrorBand(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::PlotUtils::MnvLatErrorBand*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::PlotUtils::MnvLatErrorBand

namespace PlotUtils {
//______________________________________________________________________________
void MnvH1D::Streamer(TBuffer &R__b)
{
   // Stream an object of class PlotUtils::MnvH1D.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PlotUtils::MnvH1D::Class(),this);
   } else {
      R__b.WriteClassBuffer(PlotUtils::MnvH1D::Class(),this);
   }
}

} // namespace PlotUtils
namespace ROOT {
   // Wrappers around operator new
   static void *new_PlotUtilscLcLMnvH1D(void *p) {
      return  p ? new(p) ::PlotUtils::MnvH1D : new ::PlotUtils::MnvH1D;
   }
   static void *newArray_PlotUtilscLcLMnvH1D(Long_t nElements, void *p) {
      return p ? new(p) ::PlotUtils::MnvH1D[nElements] : new ::PlotUtils::MnvH1D[nElements];
   }
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLMnvH1D(void *p) {
      delete (static_cast<::PlotUtils::MnvH1D*>(p));
   }
   static void deleteArray_PlotUtilscLcLMnvH1D(void *p) {
      delete [] (static_cast<::PlotUtils::MnvH1D*>(p));
   }
   static void destruct_PlotUtilscLcLMnvH1D(void *p) {
      typedef ::PlotUtils::MnvH1D current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
   // Wrapper around the directory auto add.
   static void directoryAutoAdd_PlotUtilscLcLMnvH1D(void *p, TDirectory *dir) {
      ((::PlotUtils::MnvH1D*)p)->DirectoryAutoAdd(dir);
   }
   // Wrapper around the merge function.
   static Long64_t  merge_PlotUtilscLcLMnvH1D(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::PlotUtils::MnvH1D*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::PlotUtils::MnvH1D

namespace ROOT {
   // Wrappers around operator new
   static void *new_PlotUtilscLcLMnvHist(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::PlotUtils::MnvHist : new ::PlotUtils::MnvHist;
   }
   static void *newArray_PlotUtilscLcLMnvHist(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::PlotUtils::MnvHist[nElements] : new ::PlotUtils::MnvHist[nElements];
   }
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLMnvHist(void *p) {
      delete (static_cast<::PlotUtils::MnvHist*>(p));
   }
   static void deleteArray_PlotUtilscLcLMnvHist(void *p) {
      delete [] (static_cast<::PlotUtils::MnvHist*>(p));
   }
   static void destruct_PlotUtilscLcLMnvHist(void *p) {
      typedef ::PlotUtils::MnvHist current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::PlotUtils::MnvHist

namespace PlotUtils {
//______________________________________________________________________________
void MnvVertErrorBand2D::Streamer(TBuffer &R__b)
{
   // Stream an object of class PlotUtils::MnvVertErrorBand2D.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PlotUtils::MnvVertErrorBand2D::Class(),this);
   } else {
      R__b.WriteClassBuffer(PlotUtils::MnvVertErrorBand2D::Class(),this);
   }
}

} // namespace PlotUtils
namespace ROOT {
   // Wrappers around operator new
   static void *new_PlotUtilscLcLMnvVertErrorBand2D(void *p) {
      return  p ? new(p) ::PlotUtils::MnvVertErrorBand2D : new ::PlotUtils::MnvVertErrorBand2D;
   }
   static void *newArray_PlotUtilscLcLMnvVertErrorBand2D(Long_t nElements, void *p) {
      return p ? new(p) ::PlotUtils::MnvVertErrorBand2D[nElements] : new ::PlotUtils::MnvVertErrorBand2D[nElements];
   }
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLMnvVertErrorBand2D(void *p) {
      delete (static_cast<::PlotUtils::MnvVertErrorBand2D*>(p));
   }
   static void deleteArray_PlotUtilscLcLMnvVertErrorBand2D(void *p) {
      delete [] (static_cast<::PlotUtils::MnvVertErrorBand2D*>(p));
   }
   static void destruct_PlotUtilscLcLMnvVertErrorBand2D(void *p) {
      typedef ::PlotUtils::MnvVertErrorBand2D current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
   // Wrapper around the directory auto add.
   static void directoryAutoAdd_PlotUtilscLcLMnvVertErrorBand2D(void *p, TDirectory *dir) {
      ((::PlotUtils::MnvVertErrorBand2D*)p)->DirectoryAutoAdd(dir);
   }
   // Wrapper around the merge function.
   static Long64_t  merge_PlotUtilscLcLMnvVertErrorBand2D(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::PlotUtils::MnvVertErrorBand2D*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::PlotUtils::MnvVertErrorBand2D

namespace PlotUtils {
//______________________________________________________________________________
void MnvLatErrorBand2D::Streamer(TBuffer &R__b)
{
   // Stream an object of class PlotUtils::MnvLatErrorBand2D.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PlotUtils::MnvLatErrorBand2D::Class(),this);
   } else {
      R__b.WriteClassBuffer(PlotUtils::MnvLatErrorBand2D::Class(),this);
   }
}

} // namespace PlotUtils
namespace ROOT {
   // Wrappers around operator new
   static void *new_PlotUtilscLcLMnvLatErrorBand2D(void *p) {
      return  p ? new(p) ::PlotUtils::MnvLatErrorBand2D : new ::PlotUtils::MnvLatErrorBand2D;
   }
   static void *newArray_PlotUtilscLcLMnvLatErrorBand2D(Long_t nElements, void *p) {
      return p ? new(p) ::PlotUtils::MnvLatErrorBand2D[nElements] : new ::PlotUtils::MnvLatErrorBand2D[nElements];
   }
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLMnvLatErrorBand2D(void *p) {
      delete (static_cast<::PlotUtils::MnvLatErrorBand2D*>(p));
   }
   static void deleteArray_PlotUtilscLcLMnvLatErrorBand2D(void *p) {
      delete [] (static_cast<::PlotUtils::MnvLatErrorBand2D*>(p));
   }
   static void destruct_PlotUtilscLcLMnvLatErrorBand2D(void *p) {
      typedef ::PlotUtils::MnvLatErrorBand2D current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
   // Wrapper around the directory auto add.
   static void directoryAutoAdd_PlotUtilscLcLMnvLatErrorBand2D(void *p, TDirectory *dir) {
      ((::PlotUtils::MnvLatErrorBand2D*)p)->DirectoryAutoAdd(dir);
   }
   // Wrapper around the merge function.
   static Long64_t  merge_PlotUtilscLcLMnvLatErrorBand2D(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::PlotUtils::MnvLatErrorBand2D*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::PlotUtils::MnvLatErrorBand2D

namespace PlotUtils {
//______________________________________________________________________________
void MnvH2D::Streamer(TBuffer &R__b)
{
   // Stream an object of class PlotUtils::MnvH2D.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PlotUtils::MnvH2D::Class(),this);
   } else {
      R__b.WriteClassBuffer(PlotUtils::MnvH2D::Class(),this);
   }
}

} // namespace PlotUtils
namespace ROOT {
   // Wrappers around operator new
   static void *new_PlotUtilscLcLMnvH2D(void *p) {
      return  p ? new(p) ::PlotUtils::MnvH2D : new ::PlotUtils::MnvH2D;
   }
   static void *newArray_PlotUtilscLcLMnvH2D(Long_t nElements, void *p) {
      return p ? new(p) ::PlotUtils::MnvH2D[nElements] : new ::PlotUtils::MnvH2D[nElements];
   }
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLMnvH2D(void *p) {
      delete (static_cast<::PlotUtils::MnvH2D*>(p));
   }
   static void deleteArray_PlotUtilscLcLMnvH2D(void *p) {
      delete [] (static_cast<::PlotUtils::MnvH2D*>(p));
   }
   static void destruct_PlotUtilscLcLMnvH2D(void *p) {
      typedef ::PlotUtils::MnvH2D current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
   // Wrapper around the directory auto add.
   static void directoryAutoAdd_PlotUtilscLcLMnvH2D(void *p, TDirectory *dir) {
      ((::PlotUtils::MnvH2D*)p)->DirectoryAutoAdd(dir);
   }
   // Wrapper around the merge function.
   static Long64_t  merge_PlotUtilscLcLMnvH2D(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::PlotUtils::MnvH2D*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::PlotUtils::MnvH2D

namespace PlotUtils {
//______________________________________________________________________________
void MnvVertErrorBand3D::Streamer(TBuffer &R__b)
{
   // Stream an object of class PlotUtils::MnvVertErrorBand3D.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PlotUtils::MnvVertErrorBand3D::Class(),this);
   } else {
      R__b.WriteClassBuffer(PlotUtils::MnvVertErrorBand3D::Class(),this);
   }
}

} // namespace PlotUtils
namespace ROOT {
   // Wrappers around operator new
   static void *new_PlotUtilscLcLMnvVertErrorBand3D(void *p) {
      return  p ? new(p) ::PlotUtils::MnvVertErrorBand3D : new ::PlotUtils::MnvVertErrorBand3D;
   }
   static void *newArray_PlotUtilscLcLMnvVertErrorBand3D(Long_t nElements, void *p) {
      return p ? new(p) ::PlotUtils::MnvVertErrorBand3D[nElements] : new ::PlotUtils::MnvVertErrorBand3D[nElements];
   }
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLMnvVertErrorBand3D(void *p) {
      delete (static_cast<::PlotUtils::MnvVertErrorBand3D*>(p));
   }
   static void deleteArray_PlotUtilscLcLMnvVertErrorBand3D(void *p) {
      delete [] (static_cast<::PlotUtils::MnvVertErrorBand3D*>(p));
   }
   static void destruct_PlotUtilscLcLMnvVertErrorBand3D(void *p) {
      typedef ::PlotUtils::MnvVertErrorBand3D current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
   // Wrapper around the directory auto add.
   static void directoryAutoAdd_PlotUtilscLcLMnvVertErrorBand3D(void *p, TDirectory *dir) {
      ((::PlotUtils::MnvVertErrorBand3D*)p)->DirectoryAutoAdd(dir);
   }
   // Wrapper around the merge function.
   static Long64_t  merge_PlotUtilscLcLMnvVertErrorBand3D(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::PlotUtils::MnvVertErrorBand3D*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::PlotUtils::MnvVertErrorBand3D

namespace PlotUtils {
//______________________________________________________________________________
void MnvLatErrorBand3D::Streamer(TBuffer &R__b)
{
   // Stream an object of class PlotUtils::MnvLatErrorBand3D.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PlotUtils::MnvLatErrorBand3D::Class(),this);
   } else {
      R__b.WriteClassBuffer(PlotUtils::MnvLatErrorBand3D::Class(),this);
   }
}

} // namespace PlotUtils
namespace ROOT {
   // Wrappers around operator new
   static void *new_PlotUtilscLcLMnvLatErrorBand3D(void *p) {
      return  p ? new(p) ::PlotUtils::MnvLatErrorBand3D : new ::PlotUtils::MnvLatErrorBand3D;
   }
   static void *newArray_PlotUtilscLcLMnvLatErrorBand3D(Long_t nElements, void *p) {
      return p ? new(p) ::PlotUtils::MnvLatErrorBand3D[nElements] : new ::PlotUtils::MnvLatErrorBand3D[nElements];
   }
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLMnvLatErrorBand3D(void *p) {
      delete (static_cast<::PlotUtils::MnvLatErrorBand3D*>(p));
   }
   static void deleteArray_PlotUtilscLcLMnvLatErrorBand3D(void *p) {
      delete [] (static_cast<::PlotUtils::MnvLatErrorBand3D*>(p));
   }
   static void destruct_PlotUtilscLcLMnvLatErrorBand3D(void *p) {
      typedef ::PlotUtils::MnvLatErrorBand3D current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
   // Wrapper around the directory auto add.
   static void directoryAutoAdd_PlotUtilscLcLMnvLatErrorBand3D(void *p, TDirectory *dir) {
      ((::PlotUtils::MnvLatErrorBand3D*)p)->DirectoryAutoAdd(dir);
   }
   // Wrapper around the merge function.
   static Long64_t  merge_PlotUtilscLcLMnvLatErrorBand3D(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::PlotUtils::MnvLatErrorBand3D*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::PlotUtils::MnvLatErrorBand3D

namespace PlotUtils {
//______________________________________________________________________________
void MnvH3D::Streamer(TBuffer &R__b)
{
   // Stream an object of class PlotUtils::MnvH3D.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PlotUtils::MnvH3D::Class(),this);
   } else {
      R__b.WriteClassBuffer(PlotUtils::MnvH3D::Class(),this);
   }
}

} // namespace PlotUtils
namespace ROOT {
   // Wrappers around operator new
   static void *new_PlotUtilscLcLMnvH3D(void *p) {
      return  p ? new(p) ::PlotUtils::MnvH3D : new ::PlotUtils::MnvH3D;
   }
   static void *newArray_PlotUtilscLcLMnvH3D(Long_t nElements, void *p) {
      return p ? new(p) ::PlotUtils::MnvH3D[nElements] : new ::PlotUtils::MnvH3D[nElements];
   }
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLMnvH3D(void *p) {
      delete (static_cast<::PlotUtils::MnvH3D*>(p));
   }
   static void deleteArray_PlotUtilscLcLMnvH3D(void *p) {
      delete [] (static_cast<::PlotUtils::MnvH3D*>(p));
   }
   static void destruct_PlotUtilscLcLMnvH3D(void *p) {
      typedef ::PlotUtils::MnvH3D current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
   // Wrapper around the directory auto add.
   static void directoryAutoAdd_PlotUtilscLcLMnvH3D(void *p, TDirectory *dir) {
      ((::PlotUtils::MnvH3D*)p)->DirectoryAutoAdd(dir);
   }
   // Wrapper around the merge function.
   static Long64_t  merge_PlotUtilscLcLMnvH3D(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::PlotUtils::MnvH3D*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::PlotUtils::MnvH3D

namespace ROOT {
   // Wrappers around operator new
   static void *new_PlotUtilscLcLMnvPlotter(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::PlotUtils::MnvPlotter : new ::PlotUtils::MnvPlotter;
   }
   static void *newArray_PlotUtilscLcLMnvPlotter(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::PlotUtils::MnvPlotter[nElements] : new ::PlotUtils::MnvPlotter[nElements];
   }
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLMnvPlotter(void *p) {
      delete (static_cast<::PlotUtils::MnvPlotter*>(p));
   }
   static void deleteArray_PlotUtilscLcLMnvPlotter(void *p) {
      delete [] (static_cast<::PlotUtils::MnvPlotter*>(p));
   }
   static void destruct_PlotUtilscLcLMnvPlotter(void *p) {
      typedef ::PlotUtils::MnvPlotter current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::PlotUtils::MnvPlotter

namespace PlotUtils {
//______________________________________________________________________________
void TreeWrapper::Streamer(TBuffer &R__b)
{
   // Stream an object of class PlotUtils::TreeWrapper.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PlotUtils::TreeWrapper::Class(),this);
   } else {
      R__b.WriteClassBuffer(PlotUtils::TreeWrapper::Class(),this);
   }
}

} // namespace PlotUtils
namespace ROOT {
   // Wrappers around operator new
   static void *new_PlotUtilscLcLTreeWrapper(void *p) {
      return  p ? new(p) ::PlotUtils::TreeWrapper : new ::PlotUtils::TreeWrapper;
   }
   static void *newArray_PlotUtilscLcLTreeWrapper(Long_t nElements, void *p) {
      return p ? new(p) ::PlotUtils::TreeWrapper[nElements] : new ::PlotUtils::TreeWrapper[nElements];
   }
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLTreeWrapper(void *p) {
      delete (static_cast<::PlotUtils::TreeWrapper*>(p));
   }
   static void deleteArray_PlotUtilscLcLTreeWrapper(void *p) {
      delete [] (static_cast<::PlotUtils::TreeWrapper*>(p));
   }
   static void destruct_PlotUtilscLcLTreeWrapper(void *p) {
      typedef ::PlotUtils::TreeWrapper current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::PlotUtils::TreeWrapper

namespace PlotUtils {
//______________________________________________________________________________
void ChainWrapper::Streamer(TBuffer &R__b)
{
   // Stream an object of class PlotUtils::ChainWrapper.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PlotUtils::ChainWrapper::Class(),this);
   } else {
      R__b.WriteClassBuffer(PlotUtils::ChainWrapper::Class(),this);
   }
}

} // namespace PlotUtils
namespace ROOT {
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLChainWrapper(void *p) {
      delete (static_cast<::PlotUtils::ChainWrapper*>(p));
   }
   static void deleteArray_PlotUtilscLcLChainWrapper(void *p) {
      delete [] (static_cast<::PlotUtils::ChainWrapper*>(p));
   }
   static void destruct_PlotUtilscLcLChainWrapper(void *p) {
      typedef ::PlotUtils::ChainWrapper current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::PlotUtils::ChainWrapper

namespace PlotUtils {
//______________________________________________________________________________
void GridCanvas::Streamer(TBuffer &R__b)
{
   // Stream an object of class PlotUtils::GridCanvas.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PlotUtils::GridCanvas::Class(),this);
   } else {
      R__b.WriteClassBuffer(PlotUtils::GridCanvas::Class(),this);
   }
}

} // namespace PlotUtils
namespace ROOT {
   // Wrappers around operator new
   static void *new_PlotUtilscLcLGridCanvas(void *p) {
      return  p ? new(p) ::PlotUtils::GridCanvas : new ::PlotUtils::GridCanvas;
   }
   static void *newArray_PlotUtilscLcLGridCanvas(Long_t nElements, void *p) {
      return p ? new(p) ::PlotUtils::GridCanvas[nElements] : new ::PlotUtils::GridCanvas[nElements];
   }
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLGridCanvas(void *p) {
      delete (static_cast<::PlotUtils::GridCanvas*>(p));
   }
   static void deleteArray_PlotUtilscLcLGridCanvas(void *p) {
      delete [] (static_cast<::PlotUtils::GridCanvas*>(p));
   }
   static void destruct_PlotUtilscLcLGridCanvas(void *p) {
      typedef ::PlotUtils::GridCanvas current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::PlotUtils::GridCanvas

namespace ROOT {
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLFluxReweighter(void *p) {
      delete (static_cast<::PlotUtils::FluxReweighter*>(p));
   }
   static void deleteArray_PlotUtilscLcLFluxReweighter(void *p) {
      delete [] (static_cast<::PlotUtils::FluxReweighter*>(p));
   }
   static void destruct_PlotUtilscLcLFluxReweighter(void *p) {
      typedef ::PlotUtils::FluxReweighter current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::PlotUtils::FluxReweighter

namespace ROOT {
   // Wrappers around operator new
   static void *new_PlotUtilscLcLBaseUniverse(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::PlotUtils::BaseUniverse : new ::PlotUtils::BaseUniverse;
   }
   static void *newArray_PlotUtilscLcLBaseUniverse(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::PlotUtils::BaseUniverse[nElements] : new ::PlotUtils::BaseUniverse[nElements];
   }
   // Wrapper around operator delete
   static void delete_PlotUtilscLcLBaseUniverse(void *p) {
      delete (static_cast<::PlotUtils::BaseUniverse*>(p));
   }
   static void deleteArray_PlotUtilscLcLBaseUniverse(void *p) {
      delete [] (static_cast<::PlotUtils::BaseUniverse*>(p));
   }
   static void destruct_PlotUtilscLcLBaseUniverse(void *p) {
      typedef ::PlotUtils::BaseUniverse current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::PlotUtils::BaseUniverse

namespace ROOT {
   // Wrapper around operator delete
   static void delete_ROOTcLcLexception(void *p) {
      delete (static_cast<::ROOT::exception*>(p));
   }
   static void deleteArray_ROOTcLcLexception(void *p) {
      delete [] (static_cast<::ROOT::exception*>(p));
   }
   static void destruct_ROOTcLcLexception(void *p) {
      typedef ::ROOT::exception current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ROOT::exception

namespace ROOT {
   // Wrapper around operator delete
   static void delete_ROOTcLcLerror(void *p) {
      delete (static_cast<::ROOT::error*>(p));
   }
   static void deleteArray_ROOTcLcLerror(void *p) {
      delete [] (static_cast<::ROOT::error*>(p));
   }
   static void destruct_ROOTcLcLerror(void *p) {
      typedef ::ROOT::error current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ROOT::error

namespace ROOT {
   // Wrapper around operator delete
   static void delete_ROOTcLcLwarning(void *p) {
      delete (static_cast<::ROOT::warning*>(p));
   }
   static void deleteArray_ROOTcLcLwarning(void *p) {
      delete [] (static_cast<::ROOT::warning*>(p));
   }
   static void destruct_ROOTcLcLwarning(void *p) {
      typedef ::ROOT::warning current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ROOT::warning

namespace ROOT {
   // Wrappers around operator new
   static void *new_ROOTcLcLdetailcLcLbeforeMain(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::ROOT::detail::beforeMain : new ::ROOT::detail::beforeMain;
   }
   static void *newArray_ROOTcLcLdetailcLcLbeforeMain(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::ROOT::detail::beforeMain[nElements] : new ::ROOT::detail::beforeMain[nElements];
   }
   // Wrapper around operator delete
   static void delete_ROOTcLcLdetailcLcLbeforeMain(void *p) {
      delete (static_cast<::ROOT::detail::beforeMain*>(p));
   }
   static void deleteArray_ROOTcLcLdetailcLcLbeforeMain(void *p) {
      delete [] (static_cast<::ROOT::detail::beforeMain*>(p));
   }
   static void destruct_ROOTcLcLdetailcLcLbeforeMain(void *p) {
      typedef ::ROOT::detail::beforeMain current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ROOT::detail::beforeMain

namespace ROOT {
   // Wrappers around operator new
   static void *new_pairlEstringcOvectorlEstringgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,vector<string> > : new pair<string,vector<string> >;
   }
   static void *newArray_pairlEstringcOvectorlEstringgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,vector<string> >[nElements] : new pair<string,vector<string> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_pairlEstringcOvectorlEstringgRsPgR(void *p) {
      delete (static_cast<pair<string,vector<string> >*>(p));
   }
   static void deleteArray_pairlEstringcOvectorlEstringgRsPgR(void *p) {
      delete [] (static_cast<pair<string,vector<string> >*>(p));
   }
   static void destruct_pairlEstringcOvectorlEstringgRsPgR(void *p) {
      typedef pair<string,vector<string> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class pair<string,vector<string> >

namespace ROOT {
   // Wrappers around operator new
   static void *new_pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,PlotUtils::MnvLatErrorBand*> : new pair<string,PlotUtils::MnvLatErrorBand*>;
   }
   static void *newArray_pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,PlotUtils::MnvLatErrorBand*>[nElements] : new pair<string,PlotUtils::MnvLatErrorBand*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(void *p) {
      delete (static_cast<pair<string,PlotUtils::MnvLatErrorBand*>*>(p));
   }
   static void deleteArray_pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(void *p) {
      delete [] (static_cast<pair<string,PlotUtils::MnvLatErrorBand*>*>(p));
   }
   static void destruct_pairlEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(void *p) {
      typedef pair<string,PlotUtils::MnvLatErrorBand*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class pair<string,PlotUtils::MnvLatErrorBand*>

namespace ROOT {
   // Wrappers around operator new
   static void *new_pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,PlotUtils::MnvVertErrorBand*> : new pair<string,PlotUtils::MnvVertErrorBand*>;
   }
   static void *newArray_pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,PlotUtils::MnvVertErrorBand*>[nElements] : new pair<string,PlotUtils::MnvVertErrorBand*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(void *p) {
      delete (static_cast<pair<string,PlotUtils::MnvVertErrorBand*>*>(p));
   }
   static void deleteArray_pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(void *p) {
      delete [] (static_cast<pair<string,PlotUtils::MnvVertErrorBand*>*>(p));
   }
   static void destruct_pairlEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(void *p) {
      typedef pair<string,PlotUtils::MnvVertErrorBand*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class pair<string,PlotUtils::MnvVertErrorBand*>

namespace ROOT {
   // Wrappers around operator new
   static void *new_pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,PlotUtils::MnvVertErrorBand3D*> : new pair<string,PlotUtils::MnvVertErrorBand3D*>;
   }
   static void *newArray_pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,PlotUtils::MnvVertErrorBand3D*>[nElements] : new pair<string,PlotUtils::MnvVertErrorBand3D*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(void *p) {
      delete (static_cast<pair<string,PlotUtils::MnvVertErrorBand3D*>*>(p));
   }
   static void deleteArray_pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(void *p) {
      delete [] (static_cast<pair<string,PlotUtils::MnvVertErrorBand3D*>*>(p));
   }
   static void destruct_pairlEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(void *p) {
      typedef pair<string,PlotUtils::MnvVertErrorBand3D*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class pair<string,PlotUtils::MnvVertErrorBand3D*>

namespace ROOT {
   // Wrappers around operator new
   static void *new_pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,PlotUtils::MnvLatErrorBand3D*> : new pair<string,PlotUtils::MnvLatErrorBand3D*>;
   }
   static void *newArray_pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,PlotUtils::MnvLatErrorBand3D*>[nElements] : new pair<string,PlotUtils::MnvLatErrorBand3D*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(void *p) {
      delete (static_cast<pair<string,PlotUtils::MnvLatErrorBand3D*>*>(p));
   }
   static void deleteArray_pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(void *p) {
      delete [] (static_cast<pair<string,PlotUtils::MnvLatErrorBand3D*>*>(p));
   }
   static void destruct_pairlEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(void *p) {
      typedef pair<string,PlotUtils::MnvLatErrorBand3D*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class pair<string,PlotUtils::MnvLatErrorBand3D*>

namespace ROOT {
   // Wrappers around operator new
   static void *new_pairlEstringcOTMatrixTlEdoublegRmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,TMatrixT<double>*> : new pair<string,TMatrixT<double>*>;
   }
   static void *newArray_pairlEstringcOTMatrixTlEdoublegRmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,TMatrixT<double>*>[nElements] : new pair<string,TMatrixT<double>*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_pairlEstringcOTMatrixTlEdoublegRmUgR(void *p) {
      delete (static_cast<pair<string,TMatrixT<double>*>*>(p));
   }
   static void deleteArray_pairlEstringcOTMatrixTlEdoublegRmUgR(void *p) {
      delete [] (static_cast<pair<string,TMatrixT<double>*>*>(p));
   }
   static void destruct_pairlEstringcOTMatrixTlEdoublegRmUgR(void *p) {
      typedef pair<string,TMatrixT<double>*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class pair<string,TMatrixT<double>*>

namespace ROOT {
   // Wrappers around operator new
   static void *new_pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,PlotUtils::MnvVertErrorBand2D*> : new pair<string,PlotUtils::MnvVertErrorBand2D*>;
   }
   static void *newArray_pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,PlotUtils::MnvVertErrorBand2D*>[nElements] : new pair<string,PlotUtils::MnvVertErrorBand2D*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(void *p) {
      delete (static_cast<pair<string,PlotUtils::MnvVertErrorBand2D*>*>(p));
   }
   static void deleteArray_pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(void *p) {
      delete [] (static_cast<pair<string,PlotUtils::MnvVertErrorBand2D*>*>(p));
   }
   static void destruct_pairlEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(void *p) {
      typedef pair<string,PlotUtils::MnvVertErrorBand2D*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class pair<string,PlotUtils::MnvVertErrorBand2D*>

namespace ROOT {
   // Wrappers around operator new
   static void *new_pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,PlotUtils::MnvLatErrorBand2D*> : new pair<string,PlotUtils::MnvLatErrorBand2D*>;
   }
   static void *newArray_pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,PlotUtils::MnvLatErrorBand2D*>[nElements] : new pair<string,PlotUtils::MnvLatErrorBand2D*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(void *p) {
      delete (static_cast<pair<string,PlotUtils::MnvLatErrorBand2D*>*>(p));
   }
   static void deleteArray_pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(void *p) {
      delete [] (static_cast<pair<string,PlotUtils::MnvLatErrorBand2D*>*>(p));
   }
   static void destruct_pairlEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(void *p) {
      typedef pair<string,PlotUtils::MnvLatErrorBand2D*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class pair<string,PlotUtils::MnvLatErrorBand2D*>

namespace ROOT {
   static TClass *vectorlEstringgR_Dictionary();
   static void vectorlEstringgR_TClassManip(TClass*);
   static void *new_vectorlEstringgR(void *p = nullptr);
   static void *newArray_vectorlEstringgR(Long_t size, void *p);
   static void delete_vectorlEstringgR(void *p);
   static void deleteArray_vectorlEstringgR(void *p);
   static void destruct_vectorlEstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<string>*)
   {
      vector<string> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<string>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<string>", -2, "vector", 423,
                  typeid(vector<string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEstringgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<string>) );
      instance.SetNew(&new_vectorlEstringgR);
      instance.SetNewArray(&newArray_vectorlEstringgR);
      instance.SetDelete(&delete_vectorlEstringgR);
      instance.SetDeleteArray(&deleteArray_vectorlEstringgR);
      instance.SetDestructor(&destruct_vectorlEstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<string> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<string>","std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<string>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<string>*>(nullptr))->GetClass();
      vectorlEstringgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string> : new vector<string>;
   }
   static void *newArray_vectorlEstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string>[nElements] : new vector<string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEstringgR(void *p) {
      delete (static_cast<vector<string>*>(p));
   }
   static void deleteArray_vectorlEstringgR(void *p) {
      delete [] (static_cast<vector<string>*>(p));
   }
   static void destruct_vectorlEstringgR(void *p) {
      typedef vector<string> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<string>

namespace ROOT {
   static TClass *vectorlETH3DmUgR_Dictionary();
   static void vectorlETH3DmUgR_TClassManip(TClass*);
   static void *new_vectorlETH3DmUgR(void *p = nullptr);
   static void *newArray_vectorlETH3DmUgR(Long_t size, void *p);
   static void delete_vectorlETH3DmUgR(void *p);
   static void deleteArray_vectorlETH3DmUgR(void *p);
   static void destruct_vectorlETH3DmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TH3D*>*)
   {
      vector<TH3D*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TH3D*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TH3D*>", -2, "vector", 423,
                  typeid(vector<TH3D*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETH3DmUgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TH3D*>) );
      instance.SetNew(&new_vectorlETH3DmUgR);
      instance.SetNewArray(&newArray_vectorlETH3DmUgR);
      instance.SetDelete(&delete_vectorlETH3DmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETH3DmUgR);
      instance.SetDestructor(&destruct_vectorlETH3DmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TH3D*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<TH3D*>","std::vector<TH3D*, std::allocator<TH3D*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TH3D*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETH3DmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TH3D*>*>(nullptr))->GetClass();
      vectorlETH3DmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETH3DmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETH3DmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TH3D*> : new vector<TH3D*>;
   }
   static void *newArray_vectorlETH3DmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TH3D*>[nElements] : new vector<TH3D*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETH3DmUgR(void *p) {
      delete (static_cast<vector<TH3D*>*>(p));
   }
   static void deleteArray_vectorlETH3DmUgR(void *p) {
      delete [] (static_cast<vector<TH3D*>*>(p));
   }
   static void destruct_vectorlETH3DmUgR(void *p) {
      typedef vector<TH3D*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TH3D*>

namespace ROOT {
   static TClass *vectorlETH2DmUgR_Dictionary();
   static void vectorlETH2DmUgR_TClassManip(TClass*);
   static void *new_vectorlETH2DmUgR(void *p = nullptr);
   static void *newArray_vectorlETH2DmUgR(Long_t size, void *p);
   static void delete_vectorlETH2DmUgR(void *p);
   static void deleteArray_vectorlETH2DmUgR(void *p);
   static void destruct_vectorlETH2DmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TH2D*>*)
   {
      vector<TH2D*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TH2D*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TH2D*>", -2, "vector", 423,
                  typeid(vector<TH2D*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETH2DmUgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TH2D*>) );
      instance.SetNew(&new_vectorlETH2DmUgR);
      instance.SetNewArray(&newArray_vectorlETH2DmUgR);
      instance.SetDelete(&delete_vectorlETH2DmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETH2DmUgR);
      instance.SetDestructor(&destruct_vectorlETH2DmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TH2D*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<TH2D*>","std::vector<TH2D*, std::allocator<TH2D*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TH2D*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETH2DmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TH2D*>*>(nullptr))->GetClass();
      vectorlETH2DmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETH2DmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETH2DmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TH2D*> : new vector<TH2D*>;
   }
   static void *newArray_vectorlETH2DmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TH2D*>[nElements] : new vector<TH2D*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETH2DmUgR(void *p) {
      delete (static_cast<vector<TH2D*>*>(p));
   }
   static void deleteArray_vectorlETH2DmUgR(void *p) {
      delete [] (static_cast<vector<TH2D*>*>(p));
   }
   static void destruct_vectorlETH2DmUgR(void *p) {
      typedef vector<TH2D*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TH2D*>

namespace ROOT {
   static TClass *vectorlETH1DmUgR_Dictionary();
   static void vectorlETH1DmUgR_TClassManip(TClass*);
   static void *new_vectorlETH1DmUgR(void *p = nullptr);
   static void *newArray_vectorlETH1DmUgR(Long_t size, void *p);
   static void delete_vectorlETH1DmUgR(void *p);
   static void deleteArray_vectorlETH1DmUgR(void *p);
   static void destruct_vectorlETH1DmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TH1D*>*)
   {
      vector<TH1D*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TH1D*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TH1D*>", -2, "vector", 423,
                  typeid(vector<TH1D*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETH1DmUgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TH1D*>) );
      instance.SetNew(&new_vectorlETH1DmUgR);
      instance.SetNewArray(&newArray_vectorlETH1DmUgR);
      instance.SetDelete(&delete_vectorlETH1DmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETH1DmUgR);
      instance.SetDestructor(&destruct_vectorlETH1DmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TH1D*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<TH1D*>","std::vector<TH1D*, std::allocator<TH1D*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TH1D*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETH1DmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TH1D*>*>(nullptr))->GetClass();
      vectorlETH1DmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETH1DmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETH1DmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TH1D*> : new vector<TH1D*>;
   }
   static void *newArray_vectorlETH1DmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TH1D*>[nElements] : new vector<TH1D*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETH1DmUgR(void *p) {
      delete (static_cast<vector<TH1D*>*>(p));
   }
   static void deleteArray_vectorlETH1DmUgR(void *p) {
      delete [] (static_cast<vector<TH1D*>*>(p));
   }
   static void destruct_vectorlETH1DmUgR(void *p) {
      typedef vector<TH1D*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TH1D*>

namespace ROOT {
   static TClass *maplEstringcOvectorlEstringgRsPgR_Dictionary();
   static void maplEstringcOvectorlEstringgRsPgR_TClassManip(TClass*);
   static void *new_maplEstringcOvectorlEstringgRsPgR(void *p = nullptr);
   static void *newArray_maplEstringcOvectorlEstringgRsPgR(Long_t size, void *p);
   static void delete_maplEstringcOvectorlEstringgRsPgR(void *p);
   static void deleteArray_maplEstringcOvectorlEstringgRsPgR(void *p);
   static void destruct_maplEstringcOvectorlEstringgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,vector<string> >*)
   {
      map<string,vector<string> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,vector<string> >));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,vector<string> >", -2, "", 71,
                  typeid(map<string,vector<string> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOvectorlEstringgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(map<string,vector<string> >) );
      instance.SetNew(&new_maplEstringcOvectorlEstringgRsPgR);
      instance.SetNewArray(&newArray_maplEstringcOvectorlEstringgRsPgR);
      instance.SetDelete(&delete_maplEstringcOvectorlEstringgRsPgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOvectorlEstringgRsPgR);
      instance.SetDestructor(&destruct_maplEstringcOvectorlEstringgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,vector<string> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,vector<string> >","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,vector<string> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOvectorlEstringgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,vector<string> >*>(nullptr))->GetClass();
      maplEstringcOvectorlEstringgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOvectorlEstringgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOvectorlEstringgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,vector<string> > : new map<string,vector<string> >;
   }
   static void *newArray_maplEstringcOvectorlEstringgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,vector<string> >[nElements] : new map<string,vector<string> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOvectorlEstringgRsPgR(void *p) {
      delete (static_cast<map<string,vector<string> >*>(p));
   }
   static void deleteArray_maplEstringcOvectorlEstringgRsPgR(void *p) {
      delete [] (static_cast<map<string,vector<string> >*>(p));
   }
   static void destruct_maplEstringcOvectorlEstringgRsPgR(void *p) {
      typedef map<string,vector<string> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,vector<string> >

namespace ROOT {
   static TClass *maplEstringcOTMatrixTlEdoublegRmUgR_Dictionary();
   static void maplEstringcOTMatrixTlEdoublegRmUgR_TClassManip(TClass*);
   static void *new_maplEstringcOTMatrixTlEdoublegRmUgR(void *p = nullptr);
   static void *newArray_maplEstringcOTMatrixTlEdoublegRmUgR(Long_t size, void *p);
   static void delete_maplEstringcOTMatrixTlEdoublegRmUgR(void *p);
   static void deleteArray_maplEstringcOTMatrixTlEdoublegRmUgR(void *p);
   static void destruct_maplEstringcOTMatrixTlEdoublegRmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,TMatrixT<double>*>*)
   {
      map<string,TMatrixT<double>*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,TMatrixT<double>*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,TMatrixT<double>*>", -2, "map", 100,
                  typeid(map<string,TMatrixT<double>*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOTMatrixTlEdoublegRmUgR_Dictionary, isa_proxy, 4,
                  sizeof(map<string,TMatrixT<double>*>) );
      instance.SetNew(&new_maplEstringcOTMatrixTlEdoublegRmUgR);
      instance.SetNewArray(&newArray_maplEstringcOTMatrixTlEdoublegRmUgR);
      instance.SetDelete(&delete_maplEstringcOTMatrixTlEdoublegRmUgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOTMatrixTlEdoublegRmUgR);
      instance.SetDestructor(&destruct_maplEstringcOTMatrixTlEdoublegRmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,TMatrixT<double>*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,TMatrixT<double>*>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, TMatrixT<double>*, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, TMatrixT<double>*> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,TMatrixT<double>*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOTMatrixTlEdoublegRmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,TMatrixT<double>*>*>(nullptr))->GetClass();
      maplEstringcOTMatrixTlEdoublegRmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOTMatrixTlEdoublegRmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOTMatrixTlEdoublegRmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,TMatrixT<double>*> : new map<string,TMatrixT<double>*>;
   }
   static void *newArray_maplEstringcOTMatrixTlEdoublegRmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,TMatrixT<double>*>[nElements] : new map<string,TMatrixT<double>*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOTMatrixTlEdoublegRmUgR(void *p) {
      delete (static_cast<map<string,TMatrixT<double>*>*>(p));
   }
   static void deleteArray_maplEstringcOTMatrixTlEdoublegRmUgR(void *p) {
      delete [] (static_cast<map<string,TMatrixT<double>*>*>(p));
   }
   static void destruct_maplEstringcOTMatrixTlEdoublegRmUgR(void *p) {
      typedef map<string,TMatrixT<double>*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,TMatrixT<double>*>

namespace ROOT {
   static TClass *maplEstringcOTH1DmUgR_Dictionary();
   static void maplEstringcOTH1DmUgR_TClassManip(TClass*);
   static void *new_maplEstringcOTH1DmUgR(void *p = nullptr);
   static void *newArray_maplEstringcOTH1DmUgR(Long_t size, void *p);
   static void delete_maplEstringcOTH1DmUgR(void *p);
   static void deleteArray_maplEstringcOTH1DmUgR(void *p);
   static void destruct_maplEstringcOTH1DmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,TH1D*>*)
   {
      map<string,TH1D*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,TH1D*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,TH1D*>", -2, "map", 100,
                  typeid(map<string,TH1D*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOTH1DmUgR_Dictionary, isa_proxy, 4,
                  sizeof(map<string,TH1D*>) );
      instance.SetNew(&new_maplEstringcOTH1DmUgR);
      instance.SetNewArray(&newArray_maplEstringcOTH1DmUgR);
      instance.SetDelete(&delete_maplEstringcOTH1DmUgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOTH1DmUgR);
      instance.SetDestructor(&destruct_maplEstringcOTH1DmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,TH1D*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,TH1D*>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, TH1D*, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, TH1D*> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,TH1D*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOTH1DmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,TH1D*>*>(nullptr))->GetClass();
      maplEstringcOTH1DmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOTH1DmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOTH1DmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,TH1D*> : new map<string,TH1D*>;
   }
   static void *newArray_maplEstringcOTH1DmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,TH1D*>[nElements] : new map<string,TH1D*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOTH1DmUgR(void *p) {
      delete (static_cast<map<string,TH1D*>*>(p));
   }
   static void deleteArray_maplEstringcOTH1DmUgR(void *p) {
      delete [] (static_cast<map<string,TH1D*>*>(p));
   }
   static void destruct_maplEstringcOTH1DmUgR(void *p) {
      typedef map<string,TH1D*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,TH1D*>

namespace ROOT {
   static TClass *maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR_Dictionary();
   static void maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR_TClassManip(TClass*);
   static void *new_maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(void *p = nullptr);
   static void *newArray_maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(Long_t size, void *p);
   static void delete_maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(void *p);
   static void deleteArray_maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(void *p);
   static void destruct_maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,PlotUtils::MnvVertErrorBand3D*>*)
   {
      map<string,PlotUtils::MnvVertErrorBand3D*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,PlotUtils::MnvVertErrorBand3D*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,PlotUtils::MnvVertErrorBand3D*>", -2, "map", 100,
                  typeid(map<string,PlotUtils::MnvVertErrorBand3D*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR_Dictionary, isa_proxy, 4,
                  sizeof(map<string,PlotUtils::MnvVertErrorBand3D*>) );
      instance.SetNew(&new_maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR);
      instance.SetNewArray(&newArray_maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR);
      instance.SetDelete(&delete_maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR);
      instance.SetDestructor(&destruct_maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,PlotUtils::MnvVertErrorBand3D*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,PlotUtils::MnvVertErrorBand3D*>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, PlotUtils::MnvVertErrorBand3D*, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, PlotUtils::MnvVertErrorBand3D*> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,PlotUtils::MnvVertErrorBand3D*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,PlotUtils::MnvVertErrorBand3D*>*>(nullptr))->GetClass();
      maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,PlotUtils::MnvVertErrorBand3D*> : new map<string,PlotUtils::MnvVertErrorBand3D*>;
   }
   static void *newArray_maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,PlotUtils::MnvVertErrorBand3D*>[nElements] : new map<string,PlotUtils::MnvVertErrorBand3D*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(void *p) {
      delete (static_cast<map<string,PlotUtils::MnvVertErrorBand3D*>*>(p));
   }
   static void deleteArray_maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(void *p) {
      delete [] (static_cast<map<string,PlotUtils::MnvVertErrorBand3D*>*>(p));
   }
   static void destruct_maplEstringcOPlotUtilscLcLMnvVertErrorBand3DmUgR(void *p) {
      typedef map<string,PlotUtils::MnvVertErrorBand3D*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,PlotUtils::MnvVertErrorBand3D*>

namespace ROOT {
   static TClass *maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR_Dictionary();
   static void maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR_TClassManip(TClass*);
   static void *new_maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(void *p = nullptr);
   static void *newArray_maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(Long_t size, void *p);
   static void delete_maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(void *p);
   static void deleteArray_maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(void *p);
   static void destruct_maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,PlotUtils::MnvVertErrorBand2D*>*)
   {
      map<string,PlotUtils::MnvVertErrorBand2D*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,PlotUtils::MnvVertErrorBand2D*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,PlotUtils::MnvVertErrorBand2D*>", -2, "map", 100,
                  typeid(map<string,PlotUtils::MnvVertErrorBand2D*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR_Dictionary, isa_proxy, 4,
                  sizeof(map<string,PlotUtils::MnvVertErrorBand2D*>) );
      instance.SetNew(&new_maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR);
      instance.SetNewArray(&newArray_maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR);
      instance.SetDelete(&delete_maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR);
      instance.SetDestructor(&destruct_maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,PlotUtils::MnvVertErrorBand2D*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,PlotUtils::MnvVertErrorBand2D*>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, PlotUtils::MnvVertErrorBand2D*, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, PlotUtils::MnvVertErrorBand2D*> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,PlotUtils::MnvVertErrorBand2D*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,PlotUtils::MnvVertErrorBand2D*>*>(nullptr))->GetClass();
      maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,PlotUtils::MnvVertErrorBand2D*> : new map<string,PlotUtils::MnvVertErrorBand2D*>;
   }
   static void *newArray_maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,PlotUtils::MnvVertErrorBand2D*>[nElements] : new map<string,PlotUtils::MnvVertErrorBand2D*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(void *p) {
      delete (static_cast<map<string,PlotUtils::MnvVertErrorBand2D*>*>(p));
   }
   static void deleteArray_maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(void *p) {
      delete [] (static_cast<map<string,PlotUtils::MnvVertErrorBand2D*>*>(p));
   }
   static void destruct_maplEstringcOPlotUtilscLcLMnvVertErrorBand2DmUgR(void *p) {
      typedef map<string,PlotUtils::MnvVertErrorBand2D*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,PlotUtils::MnvVertErrorBand2D*>

namespace ROOT {
   static TClass *maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR_Dictionary();
   static void maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR_TClassManip(TClass*);
   static void *new_maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(void *p = nullptr);
   static void *newArray_maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(Long_t size, void *p);
   static void delete_maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(void *p);
   static void deleteArray_maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(void *p);
   static void destruct_maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,PlotUtils::MnvVertErrorBand*>*)
   {
      map<string,PlotUtils::MnvVertErrorBand*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,PlotUtils::MnvVertErrorBand*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,PlotUtils::MnvVertErrorBand*>", -2, "map", 100,
                  typeid(map<string,PlotUtils::MnvVertErrorBand*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR_Dictionary, isa_proxy, 4,
                  sizeof(map<string,PlotUtils::MnvVertErrorBand*>) );
      instance.SetNew(&new_maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR);
      instance.SetNewArray(&newArray_maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR);
      instance.SetDelete(&delete_maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR);
      instance.SetDestructor(&destruct_maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,PlotUtils::MnvVertErrorBand*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,PlotUtils::MnvVertErrorBand*>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, PlotUtils::MnvVertErrorBand*, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, PlotUtils::MnvVertErrorBand*> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,PlotUtils::MnvVertErrorBand*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,PlotUtils::MnvVertErrorBand*>*>(nullptr))->GetClass();
      maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,PlotUtils::MnvVertErrorBand*> : new map<string,PlotUtils::MnvVertErrorBand*>;
   }
   static void *newArray_maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,PlotUtils::MnvVertErrorBand*>[nElements] : new map<string,PlotUtils::MnvVertErrorBand*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(void *p) {
      delete (static_cast<map<string,PlotUtils::MnvVertErrorBand*>*>(p));
   }
   static void deleteArray_maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(void *p) {
      delete [] (static_cast<map<string,PlotUtils::MnvVertErrorBand*>*>(p));
   }
   static void destruct_maplEstringcOPlotUtilscLcLMnvVertErrorBandmUgR(void *p) {
      typedef map<string,PlotUtils::MnvVertErrorBand*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,PlotUtils::MnvVertErrorBand*>

namespace ROOT {
   static TClass *maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR_Dictionary();
   static void maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR_TClassManip(TClass*);
   static void *new_maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(void *p = nullptr);
   static void *newArray_maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(Long_t size, void *p);
   static void delete_maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(void *p);
   static void deleteArray_maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(void *p);
   static void destruct_maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,PlotUtils::MnvLatErrorBand3D*>*)
   {
      map<string,PlotUtils::MnvLatErrorBand3D*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,PlotUtils::MnvLatErrorBand3D*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,PlotUtils::MnvLatErrorBand3D*>", -2, "map", 100,
                  typeid(map<string,PlotUtils::MnvLatErrorBand3D*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR_Dictionary, isa_proxy, 4,
                  sizeof(map<string,PlotUtils::MnvLatErrorBand3D*>) );
      instance.SetNew(&new_maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR);
      instance.SetNewArray(&newArray_maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR);
      instance.SetDelete(&delete_maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR);
      instance.SetDestructor(&destruct_maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,PlotUtils::MnvLatErrorBand3D*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,PlotUtils::MnvLatErrorBand3D*>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, PlotUtils::MnvLatErrorBand3D*, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, PlotUtils::MnvLatErrorBand3D*> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,PlotUtils::MnvLatErrorBand3D*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,PlotUtils::MnvLatErrorBand3D*>*>(nullptr))->GetClass();
      maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,PlotUtils::MnvLatErrorBand3D*> : new map<string,PlotUtils::MnvLatErrorBand3D*>;
   }
   static void *newArray_maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,PlotUtils::MnvLatErrorBand3D*>[nElements] : new map<string,PlotUtils::MnvLatErrorBand3D*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(void *p) {
      delete (static_cast<map<string,PlotUtils::MnvLatErrorBand3D*>*>(p));
   }
   static void deleteArray_maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(void *p) {
      delete [] (static_cast<map<string,PlotUtils::MnvLatErrorBand3D*>*>(p));
   }
   static void destruct_maplEstringcOPlotUtilscLcLMnvLatErrorBand3DmUgR(void *p) {
      typedef map<string,PlotUtils::MnvLatErrorBand3D*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,PlotUtils::MnvLatErrorBand3D*>

namespace ROOT {
   static TClass *maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR_Dictionary();
   static void maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR_TClassManip(TClass*);
   static void *new_maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(void *p = nullptr);
   static void *newArray_maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(Long_t size, void *p);
   static void delete_maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(void *p);
   static void deleteArray_maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(void *p);
   static void destruct_maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,PlotUtils::MnvLatErrorBand2D*>*)
   {
      map<string,PlotUtils::MnvLatErrorBand2D*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,PlotUtils::MnvLatErrorBand2D*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,PlotUtils::MnvLatErrorBand2D*>", -2, "map", 100,
                  typeid(map<string,PlotUtils::MnvLatErrorBand2D*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR_Dictionary, isa_proxy, 4,
                  sizeof(map<string,PlotUtils::MnvLatErrorBand2D*>) );
      instance.SetNew(&new_maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR);
      instance.SetNewArray(&newArray_maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR);
      instance.SetDelete(&delete_maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR);
      instance.SetDestructor(&destruct_maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,PlotUtils::MnvLatErrorBand2D*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,PlotUtils::MnvLatErrorBand2D*>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, PlotUtils::MnvLatErrorBand2D*, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, PlotUtils::MnvLatErrorBand2D*> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,PlotUtils::MnvLatErrorBand2D*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,PlotUtils::MnvLatErrorBand2D*>*>(nullptr))->GetClass();
      maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,PlotUtils::MnvLatErrorBand2D*> : new map<string,PlotUtils::MnvLatErrorBand2D*>;
   }
   static void *newArray_maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,PlotUtils::MnvLatErrorBand2D*>[nElements] : new map<string,PlotUtils::MnvLatErrorBand2D*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(void *p) {
      delete (static_cast<map<string,PlotUtils::MnvLatErrorBand2D*>*>(p));
   }
   static void deleteArray_maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(void *p) {
      delete [] (static_cast<map<string,PlotUtils::MnvLatErrorBand2D*>*>(p));
   }
   static void destruct_maplEstringcOPlotUtilscLcLMnvLatErrorBand2DmUgR(void *p) {
      typedef map<string,PlotUtils::MnvLatErrorBand2D*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,PlotUtils::MnvLatErrorBand2D*>

namespace ROOT {
   static TClass *maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR_Dictionary();
   static void maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR_TClassManip(TClass*);
   static void *new_maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(void *p = nullptr);
   static void *newArray_maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(Long_t size, void *p);
   static void delete_maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(void *p);
   static void deleteArray_maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(void *p);
   static void destruct_maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,PlotUtils::MnvLatErrorBand*>*)
   {
      map<string,PlotUtils::MnvLatErrorBand*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,PlotUtils::MnvLatErrorBand*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,PlotUtils::MnvLatErrorBand*>", -2, "map", 100,
                  typeid(map<string,PlotUtils::MnvLatErrorBand*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR_Dictionary, isa_proxy, 4,
                  sizeof(map<string,PlotUtils::MnvLatErrorBand*>) );
      instance.SetNew(&new_maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR);
      instance.SetNewArray(&newArray_maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR);
      instance.SetDelete(&delete_maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR);
      instance.SetDestructor(&destruct_maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,PlotUtils::MnvLatErrorBand*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,PlotUtils::MnvLatErrorBand*>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, PlotUtils::MnvLatErrorBand*, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, PlotUtils::MnvLatErrorBand*> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,PlotUtils::MnvLatErrorBand*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,PlotUtils::MnvLatErrorBand*>*>(nullptr))->GetClass();
      maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,PlotUtils::MnvLatErrorBand*> : new map<string,PlotUtils::MnvLatErrorBand*>;
   }
   static void *newArray_maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,PlotUtils::MnvLatErrorBand*>[nElements] : new map<string,PlotUtils::MnvLatErrorBand*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(void *p) {
      delete (static_cast<map<string,PlotUtils::MnvLatErrorBand*>*>(p));
   }
   static void deleteArray_maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(void *p) {
      delete [] (static_cast<map<string,PlotUtils::MnvLatErrorBand*>*>(p));
   }
   static void destruct_maplEstringcOPlotUtilscLcLMnvLatErrorBandmUgR(void *p) {
      typedef map<string,PlotUtils::MnvLatErrorBand*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,PlotUtils::MnvLatErrorBand*>

namespace {
  void TriggerDictionaryInitialization_MAT_Impl() {
    static const char* headers[] = {
"0",
nullptr
    };
    static const char* includePaths[] = {
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "MAT dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$TH1.h")))  __attribute__((annotate("$clingAutoload$../PlotUtils/HistogramUtils.h")))  TH1D;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$TH2.h")))  __attribute__((annotate("$clingAutoload$../PlotUtils/HistogramUtils.h")))  TH2D;
class __attribute__((annotate("$clingAutoload$TH3.h")))  __attribute__((annotate("$clingAutoload$../PlotUtils/HistogramUtils.h")))  TH3D;
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$../PlotUtils/Exceptions.h")))  BadUnivWgtError;}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$../PlotUtils/Exceptions.h")))  ConstraintLoadError;}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$../PlotUtils/Exceptions.h")))  ConstraintAccessError;}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$../PlotUtils/Exceptions.h")))  MissingFluxUnivError;}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$../PlotUtils/Exceptions.h")))  MissingSpectatorStrategyError;}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$../PlotUtils/Exceptions.h")))  NoFluxUnivError;}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$../PlotUtils/Exceptions.h")))  NoWgtdSpreadError;}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$../PlotUtils/Exceptions.h")))  SpectatorConstraintCollisionError;}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$PlotUtils/MnvVertErrorBand.h")))  __attribute__((annotate("$clingAutoload$../PlotUtils/HistogramUtils.h")))  MnvVertErrorBand;}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$PlotUtils/MnvLatErrorBand.h")))  __attribute__((annotate("$clingAutoload$../PlotUtils/HistogramUtils.h")))  MnvLatErrorBand;}
namespace std{template <typename _CharT> struct __attribute__((annotate("$clingAutoload$bits/char_traits.h")))  __attribute__((annotate("$clingAutoload$string")))  char_traits;
}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$PlotUtils/MnvH1D.h")))  __attribute__((annotate("$clingAutoload$../PlotUtils/HistogramUtils.h")))  MnvH1D;}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$../PlotUtils/HistogramUtils.h")))  MnvHist;}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$../PlotUtils/MnvH2D.h")))  MnvH2D;}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$../PlotUtils/MnvH3D.h")))  MnvH3D;}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$../PlotUtils/MnvPlotter.h")))  MnvPlotter;}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$PlotUtils/TreeWrapper.h")))  __attribute__((annotate("$clingAutoload$../PlotUtils/ChainWrapper.h")))  TreeWrapper;}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$../PlotUtils/ChainWrapper.h")))  ChainWrapper;}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$../PlotUtils/GridCanvas.h")))  GridCanvas;}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$../PlotUtils/FluxReweighter.h")))  FluxReweighter;}
namespace PlotUtils{class __attribute__((annotate("$clingAutoload$../PlotUtils/BaseUniverse.h")))  BaseUniverse;}
namespace ROOT{class __attribute__((annotate("$clingAutoload$../PlotUtils/ErrorHandler.h")))  exception;}
namespace ROOT{class __attribute__((annotate("$clingAutoload$../PlotUtils/ErrorHandler.h")))  error;}
namespace ROOT{class __attribute__((annotate("$clingAutoload$../PlotUtils/ErrorHandler.h")))  warning;}
namespace ROOT{namespace detail{struct __attribute__((annotate("$clingAutoload$../PlotUtils/ErrorHandler.h")))  beforeMain;}}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "MAT dictionary payload"

#ifndef __GCC_XML__
  #define __GCC_XML__ 1
#endif
#ifndef NCINTEX
  #define NCINTEX 1
#endif
#ifndef FORM
  #define FORM 1
#endif
#ifndef MNVROOT6
  #define MNVROOT6 1
#endif
#ifndef PLOTUTILS_STANDALONE
  #define PLOTUTILS_STANDALONE 1
#endif
#ifndef BUILD_SHARED_LIBS
  #define BUILD_SHARED_LIBS OFF
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#ifndef PLOTUTILSDICT_H 
#define PLOTUTILSDICT_H 1

// Include files for PlotUtils dictionary.

/** @file PlotUtilsDict.h
 *  
 *
 *  @author Jeremy Wolcott <jwolcott@fnal.gov>
 *  @date   2012-11-25
 */
// ============================================================================
// PlotUtils
// ============================================================================

// here we need to include all the header files
// for the classes we want to make dictionaries of

#include <vector>

#include "../PlotUtils/Exceptions.h"
#include "../PlotUtils/HistogramUtils.h"
#include "../PlotUtils/MnvH1D.h"
#include "../PlotUtils/MnvH2D.h"
#include "../PlotUtils/MnvH3D.h"
#include "../PlotUtils/MnvLatErrorBand.h"
#include "../PlotUtils/MnvLatErrorBand2D.h"
#include "../PlotUtils/MnvLatErrorBand3D.h"
#include "../PlotUtils/MnvPlotter.h"
#include "../PlotUtils/MnvColors.h"
#include "../PlotUtils/MnvVertErrorBand.h"
#include "../PlotUtils/MnvVertErrorBand2D.h"
#include "../PlotUtils/MnvVertErrorBand3D.h"
#include "../PlotUtils/ChainWrapper.h"
#include "../PlotUtils/TreeWrapper.h"
#include "../PlotUtils/GridCanvas.h"
#include "../PlotUtils/FluxReweighter.h"
#include "../PlotUtils/BaseUniverse.h"

//TODO: Do I need this?
#include "../PlotUtils/ErrorHandler.h"

// this garbage is necessary so that gccxml is able to create dictionaries for these custom containers
// (since it otherwise doesn't know which specific version of these templated classes to instantiate)
// see: http://root.cern.ch/root/roottalk/roottalk10/0035.html
// somehow std::map<>s seem to be instantiated somewhere else, so explicit instantiation is not necessary?
//#ifdef __GCCXML__
//template class std::vector<PlotUtils::MnvEVD::Event>;                                       // the 'Events' typedef
//template class std::pair<std::string, std::vector<PlotUtils::MnvEVD::Event> >;              // the 'EventGroup' typedef

template class std::map<std::string, std::vector<std::string> >;
template class std::pair<std::string, std::vector<std::string> >;

// The std::pair<>s for those std::map<>s don't seem to be generated though.
template class std::pair< std::string, PlotUtils::MnvLatErrorBand* >;
template class std::pair< std::string, PlotUtils::MnvVertErrorBand* >;
template class std::pair< std::string, PlotUtils::MnvVertErrorBand3D* >;
template class std::pair< std::string, PlotUtils::MnvLatErrorBand3D* >;
template class std::pair< std::string, TMatrixT<double>* >;
template class std::pair< std::string, PlotUtils::MnvVertErrorBand2D* >;
template class std::pair< std::string, PlotUtils::MnvLatErrorBand2D* >;

//#endif
#endif // PLOTUTILSDICT_H


#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"MnvColors::GetColors", payloadCode, "@",
"PlotUtils::BadUnivWgtError", payloadCode, "@",
"PlotUtils::BaseUniverse", payloadCode, "@",
"PlotUtils::ChainWrapper", payloadCode, "@",
"PlotUtils::ConstraintAccessError", payloadCode, "@",
"PlotUtils::ConstraintLoadError", payloadCode, "@",
"PlotUtils::FluxReweighter", payloadCode, "@",
"PlotUtils::GridCanvas", payloadCode, "@",
"PlotUtils::HandleErrorsInPython", payloadCode, "@",
"PlotUtils::MissingFluxUnivError", payloadCode, "@",
"PlotUtils::MissingSpectatorStrategyError", payloadCode, "@",
"PlotUtils::MnvH1D", payloadCode, "@",
"PlotUtils::MnvH2D", payloadCode, "@",
"PlotUtils::MnvH3D", payloadCode, "@",
"PlotUtils::MnvHist", payloadCode, "@",
"PlotUtils::MnvLatErrorBand", payloadCode, "@",
"PlotUtils::MnvLatErrorBand2D", payloadCode, "@",
"PlotUtils::MnvLatErrorBand3D", payloadCode, "@",
"PlotUtils::MnvPlotter", payloadCode, "@",
"PlotUtils::MnvVertErrorBand", payloadCode, "@",
"PlotUtils::MnvVertErrorBand2D", payloadCode, "@",
"PlotUtils::MnvVertErrorBand3D", payloadCode, "@",
"PlotUtils::NoFluxUnivError", payloadCode, "@",
"PlotUtils::NoWgtdSpreadError", payloadCode, "@",
"PlotUtils::SpectatorConstraintCollisionError", payloadCode, "@",
"PlotUtils::TreeWrapper", payloadCode, "@",
"PlotUtils::flux_reweighter", payloadCode, "@",
"ROOT::detail::beforeMain", payloadCode, "@",
"ROOT::error", payloadCode, "@",
"ROOT::exception", payloadCode, "@",
"ROOT::warning", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("MAT",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_MAT_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_MAT_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_MAT() {
  TriggerDictionaryInitialization_MAT_Impl();
}
