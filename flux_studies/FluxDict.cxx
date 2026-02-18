// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME FluxDict
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
#include "FluxCalculatorLoop.h"
#include "FluxCalculatorLoop.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *FluxCalculatorLoop_Dictionary();
   static void FluxCalculatorLoop_TClassManip(TClass*);
   static void *new_FluxCalculatorLoop(void *p = nullptr);
   static void *newArray_FluxCalculatorLoop(Long_t size, void *p);
   static void delete_FluxCalculatorLoop(void *p);
   static void deleteArray_FluxCalculatorLoop(void *p);
   static void destruct_FluxCalculatorLoop(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FluxCalculatorLoop*)
   {
      ::FluxCalculatorLoop *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::FluxCalculatorLoop));
      static ::ROOT::TGenericClassInfo 
         instance("FluxCalculatorLoop", "FluxCalculatorLoop.h", 10,
                  typeid(::FluxCalculatorLoop), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &FluxCalculatorLoop_Dictionary, isa_proxy, 4,
                  sizeof(::FluxCalculatorLoop) );
      instance.SetNew(&new_FluxCalculatorLoop);
      instance.SetNewArray(&newArray_FluxCalculatorLoop);
      instance.SetDelete(&delete_FluxCalculatorLoop);
      instance.SetDeleteArray(&deleteArray_FluxCalculatorLoop);
      instance.SetDestructor(&destruct_FluxCalculatorLoop);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FluxCalculatorLoop*)
   {
      return GenerateInitInstanceLocal(static_cast<::FluxCalculatorLoop*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FluxCalculatorLoop*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *FluxCalculatorLoop_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::FluxCalculatorLoop*>(nullptr))->GetClass();
      FluxCalculatorLoop_TClassManip(theClass);
   return theClass;
   }

   static void FluxCalculatorLoop_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_FluxCalculatorLoop(void *p) {
      return  p ? new(p) ::FluxCalculatorLoop : new ::FluxCalculatorLoop;
   }
   static void *newArray_FluxCalculatorLoop(Long_t nElements, void *p) {
      return p ? new(p) ::FluxCalculatorLoop[nElements] : new ::FluxCalculatorLoop[nElements];
   }
   // Wrapper around operator delete
   static void delete_FluxCalculatorLoop(void *p) {
      delete (static_cast<::FluxCalculatorLoop*>(p));
   }
   static void deleteArray_FluxCalculatorLoop(void *p) {
      delete [] (static_cast<::FluxCalculatorLoop*>(p));
   }
   static void destruct_FluxCalculatorLoop(void *p) {
      typedef ::FluxCalculatorLoop current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FluxCalculatorLoop

namespace {
  void TriggerDictionaryInitialization_FluxDict_Impl() {
    static const char* headers[] = {
"FluxCalculatorLoop.h",
nullptr
    };
    static const char* includePaths[] = {
"////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/opt/spack/linux-almalinux9-x86_64_v3/gcc-12.2.0/root-6.28.12-hljl7gyomotoqht2uzvhnf73337jq67q/include/root",
"/exp/minerva/app/users/qvuong/MAT_AL9/CC-NuE-XSec/flux_studies/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "FluxDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$FluxCalculatorLoop.h")))  FluxCalculatorLoop;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "FluxDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "FluxCalculatorLoop.h"
// LinkDef.h
#pragma once
#include "FluxCalculatorLoop.h"

#ifdef __CINT__
#pragma link C++ class FluxCalculatorLoop+;
#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"FluxCalculatorLoop", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("FluxDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_FluxDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_FluxDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_FluxDict() {
  TriggerDictionaryInitialization_FluxDict_Impl();
}
