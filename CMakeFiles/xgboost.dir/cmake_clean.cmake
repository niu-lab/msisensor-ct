file(REMOVE_RECURSE
  "vendor/xgboost/lib/libxgboost.pdb"
  "vendor/xgboost/lib/libxgboost.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/xgboost.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
