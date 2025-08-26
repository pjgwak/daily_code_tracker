// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME objdIdict_Analysis
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
#include "src/Final2DFit.h"
#include "src/CtauBkgFit.h"
#include "src/CtauTrueFit.h"
#include "src/CtauResFit.h"
#include "src/CtauErrFit.h"
#include "src/MassFit.h"
#include "src/McMassFit.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *Final2DFit_Dictionary();
   static void Final2DFit_TClassManip(TClass*);
   static void delete_Final2DFit(void *p);
   static void deleteArray_Final2DFit(void *p);
   static void destruct_Final2DFit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Final2DFit*)
   {
      ::Final2DFit *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Final2DFit));
      static ::ROOT::TGenericClassInfo 
         instance("Final2DFit", "Final2DFit.h", 8,
                  typeid(::Final2DFit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Final2DFit_Dictionary, isa_proxy, 4,
                  sizeof(::Final2DFit) );
      instance.SetDelete(&delete_Final2DFit);
      instance.SetDeleteArray(&deleteArray_Final2DFit);
      instance.SetDestructor(&destruct_Final2DFit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Final2DFit*)
   {
      return GenerateInitInstanceLocal(static_cast<::Final2DFit*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Final2DFit*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Final2DFit_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Final2DFit*>(nullptr))->GetClass();
      Final2DFit_TClassManip(theClass);
   return theClass;
   }

   static void Final2DFit_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *CtauBkgFit_Dictionary();
   static void CtauBkgFit_TClassManip(TClass*);
   static void delete_CtauBkgFit(void *p);
   static void deleteArray_CtauBkgFit(void *p);
   static void destruct_CtauBkgFit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CtauBkgFit*)
   {
      ::CtauBkgFit *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::CtauBkgFit));
      static ::ROOT::TGenericClassInfo 
         instance("CtauBkgFit", "CtauBkgFit.h", 8,
                  typeid(::CtauBkgFit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &CtauBkgFit_Dictionary, isa_proxy, 4,
                  sizeof(::CtauBkgFit) );
      instance.SetDelete(&delete_CtauBkgFit);
      instance.SetDeleteArray(&deleteArray_CtauBkgFit);
      instance.SetDestructor(&destruct_CtauBkgFit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CtauBkgFit*)
   {
      return GenerateInitInstanceLocal(static_cast<::CtauBkgFit*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::CtauBkgFit*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CtauBkgFit_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::CtauBkgFit*>(nullptr))->GetClass();
      CtauBkgFit_TClassManip(theClass);
   return theClass;
   }

   static void CtauBkgFit_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *CtauTrueFit_Dictionary();
   static void CtauTrueFit_TClassManip(TClass*);
   static void delete_CtauTrueFit(void *p);
   static void deleteArray_CtauTrueFit(void *p);
   static void destruct_CtauTrueFit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CtauTrueFit*)
   {
      ::CtauTrueFit *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::CtauTrueFit));
      static ::ROOT::TGenericClassInfo 
         instance("CtauTrueFit", "CtauTrueFit.h", 8,
                  typeid(::CtauTrueFit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &CtauTrueFit_Dictionary, isa_proxy, 4,
                  sizeof(::CtauTrueFit) );
      instance.SetDelete(&delete_CtauTrueFit);
      instance.SetDeleteArray(&deleteArray_CtauTrueFit);
      instance.SetDestructor(&destruct_CtauTrueFit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CtauTrueFit*)
   {
      return GenerateInitInstanceLocal(static_cast<::CtauTrueFit*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::CtauTrueFit*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CtauTrueFit_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::CtauTrueFit*>(nullptr))->GetClass();
      CtauTrueFit_TClassManip(theClass);
   return theClass;
   }

   static void CtauTrueFit_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *CtauResFit_Dictionary();
   static void CtauResFit_TClassManip(TClass*);
   static void delete_CtauResFit(void *p);
   static void deleteArray_CtauResFit(void *p);
   static void destruct_CtauResFit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CtauResFit*)
   {
      ::CtauResFit *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::CtauResFit));
      static ::ROOT::TGenericClassInfo 
         instance("CtauResFit", "CtauResFit.h", 8,
                  typeid(::CtauResFit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &CtauResFit_Dictionary, isa_proxy, 4,
                  sizeof(::CtauResFit) );
      instance.SetDelete(&delete_CtauResFit);
      instance.SetDeleteArray(&deleteArray_CtauResFit);
      instance.SetDestructor(&destruct_CtauResFit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CtauResFit*)
   {
      return GenerateInitInstanceLocal(static_cast<::CtauResFit*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::CtauResFit*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CtauResFit_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::CtauResFit*>(nullptr))->GetClass();
      CtauResFit_TClassManip(theClass);
   return theClass;
   }

   static void CtauResFit_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *CtauErrFit_Dictionary();
   static void CtauErrFit_TClassManip(TClass*);
   static void delete_CtauErrFit(void *p);
   static void deleteArray_CtauErrFit(void *p);
   static void destruct_CtauErrFit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CtauErrFit*)
   {
      ::CtauErrFit *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::CtauErrFit));
      static ::ROOT::TGenericClassInfo 
         instance("CtauErrFit", "CtauErrFit.h", 8,
                  typeid(::CtauErrFit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &CtauErrFit_Dictionary, isa_proxy, 4,
                  sizeof(::CtauErrFit) );
      instance.SetDelete(&delete_CtauErrFit);
      instance.SetDeleteArray(&deleteArray_CtauErrFit);
      instance.SetDestructor(&destruct_CtauErrFit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CtauErrFit*)
   {
      return GenerateInitInstanceLocal(static_cast<::CtauErrFit*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::CtauErrFit*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CtauErrFit_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::CtauErrFit*>(nullptr))->GetClass();
      CtauErrFit_TClassManip(theClass);
   return theClass;
   }

   static void CtauErrFit_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *MassFit_Dictionary();
   static void MassFit_TClassManip(TClass*);
   static void delete_MassFit(void *p);
   static void deleteArray_MassFit(void *p);
   static void destruct_MassFit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MassFit*)
   {
      ::MassFit *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::MassFit));
      static ::ROOT::TGenericClassInfo 
         instance("MassFit", "MassFit.h", 8,
                  typeid(::MassFit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &MassFit_Dictionary, isa_proxy, 4,
                  sizeof(::MassFit) );
      instance.SetDelete(&delete_MassFit);
      instance.SetDeleteArray(&deleteArray_MassFit);
      instance.SetDestructor(&destruct_MassFit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MassFit*)
   {
      return GenerateInitInstanceLocal(static_cast<::MassFit*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::MassFit*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *MassFit_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::MassFit*>(nullptr))->GetClass();
      MassFit_TClassManip(theClass);
   return theClass;
   }

   static void MassFit_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *McMassFit_Dictionary();
   static void McMassFit_TClassManip(TClass*);
   static void delete_McMassFit(void *p);
   static void deleteArray_McMassFit(void *p);
   static void destruct_McMassFit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::McMassFit*)
   {
      ::McMassFit *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::McMassFit));
      static ::ROOT::TGenericClassInfo 
         instance("McMassFit", "McMassFit.h", 8,
                  typeid(::McMassFit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &McMassFit_Dictionary, isa_proxy, 4,
                  sizeof(::McMassFit) );
      instance.SetDelete(&delete_McMassFit);
      instance.SetDeleteArray(&deleteArray_McMassFit);
      instance.SetDestructor(&destruct_McMassFit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::McMassFit*)
   {
      return GenerateInitInstanceLocal(static_cast<::McMassFit*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::McMassFit*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *McMassFit_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::McMassFit*>(nullptr))->GetClass();
      McMassFit_TClassManip(theClass);
   return theClass;
   }

   static void McMassFit_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Final2DFit(void *p) {
      delete (static_cast<::Final2DFit*>(p));
   }
   static void deleteArray_Final2DFit(void *p) {
      delete [] (static_cast<::Final2DFit*>(p));
   }
   static void destruct_Final2DFit(void *p) {
      typedef ::Final2DFit current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Final2DFit

namespace ROOT {
   // Wrapper around operator delete
   static void delete_CtauBkgFit(void *p) {
      delete (static_cast<::CtauBkgFit*>(p));
   }
   static void deleteArray_CtauBkgFit(void *p) {
      delete [] (static_cast<::CtauBkgFit*>(p));
   }
   static void destruct_CtauBkgFit(void *p) {
      typedef ::CtauBkgFit current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::CtauBkgFit

namespace ROOT {
   // Wrapper around operator delete
   static void delete_CtauTrueFit(void *p) {
      delete (static_cast<::CtauTrueFit*>(p));
   }
   static void deleteArray_CtauTrueFit(void *p) {
      delete [] (static_cast<::CtauTrueFit*>(p));
   }
   static void destruct_CtauTrueFit(void *p) {
      typedef ::CtauTrueFit current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::CtauTrueFit

namespace ROOT {
   // Wrapper around operator delete
   static void delete_CtauResFit(void *p) {
      delete (static_cast<::CtauResFit*>(p));
   }
   static void deleteArray_CtauResFit(void *p) {
      delete [] (static_cast<::CtauResFit*>(p));
   }
   static void destruct_CtauResFit(void *p) {
      typedef ::CtauResFit current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::CtauResFit

namespace ROOT {
   // Wrapper around operator delete
   static void delete_CtauErrFit(void *p) {
      delete (static_cast<::CtauErrFit*>(p));
   }
   static void deleteArray_CtauErrFit(void *p) {
      delete [] (static_cast<::CtauErrFit*>(p));
   }
   static void destruct_CtauErrFit(void *p) {
      typedef ::CtauErrFit current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::CtauErrFit

namespace ROOT {
   // Wrapper around operator delete
   static void delete_MassFit(void *p) {
      delete (static_cast<::MassFit*>(p));
   }
   static void deleteArray_MassFit(void *p) {
      delete [] (static_cast<::MassFit*>(p));
   }
   static void destruct_MassFit(void *p) {
      typedef ::MassFit current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::MassFit

namespace ROOT {
   // Wrapper around operator delete
   static void delete_McMassFit(void *p) {
      delete (static_cast<::McMassFit*>(p));
   }
   static void deleteArray_McMassFit(void *p) {
      delete [] (static_cast<::McMassFit*>(p));
   }
   static void destruct_McMassFit(void *p) {
      typedef ::McMassFit current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::McMassFit

namespace {
  void TriggerDictionaryInitialization_dict_Analysis_Impl() {
    static const char* headers[] = {
"src/Final2DFit.h",
"src/CtauBkgFit.h",
"src/CtauTrueFit.h",
"src/CtauResFit.h",
"src/CtauErrFit.h",
"src/MassFit.h",
"src/McMassFit.h",
nullptr
    };
    static const char* includePaths[] = {
"src",
"/home/hep319/anaconda3/envs/root632/include/",
"/work/pjgwak/daily_code_tracker/2025/250813/run3_pbpb_mass_counting/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "dict_Analysis dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$src/Final2DFit.h")))  Final2DFit;
class __attribute__((annotate("$clingAutoload$src/CtauBkgFit.h")))  CtauBkgFit;
class __attribute__((annotate("$clingAutoload$src/CtauTrueFit.h")))  CtauTrueFit;
class __attribute__((annotate("$clingAutoload$src/CtauResFit.h")))  CtauResFit;
class __attribute__((annotate("$clingAutoload$src/CtauErrFit.h")))  CtauErrFit;
class __attribute__((annotate("$clingAutoload$src/MassFit.h")))  MassFit;
class __attribute__((annotate("$clingAutoload$src/McMassFit.h")))  McMassFit;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "dict_Analysis dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "src/Final2DFit.h"
#include "src/CtauBkgFit.h"
#include "src/CtauTrueFit.h"
#include "src/CtauResFit.h"
#include "src/CtauErrFit.h"
#include "src/MassFit.h"
#include "src/McMassFit.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"CtauBkgFit", payloadCode, "@",
"CtauErrFit", payloadCode, "@",
"CtauResFit", payloadCode, "@",
"CtauTrueFit", payloadCode, "@",
"Final2DFit", payloadCode, "@",
"MassFit", payloadCode, "@",
"McMassFit", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("dict_Analysis",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_dict_Analysis_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_dict_Analysis_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_dict_Analysis() {
  TriggerDictionaryInitialization_dict_Analysis_Impl();
}
