#ifdef _MSC_VER
  #ifdef BUILD_DLL
    #define API __declspec(dllexport)
  #else
    #define API // __declspec(dllimport)
  #endif
#else
  #define API 
#endif
