;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:1 within `advance_check2`
define void @julia_advance_check2_2812({}* nonnull align 8 dereferenceable(24) %0, { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }* nocapture nonnull readonly align 8 dereferenceable(344) %1, { double, double, double, double, double, double, double, double, {}*, {}* }* nocapture nonnull readonly align 8 dereferenceable(80) %2, { i64, i64, i64, i64, i64, i64, i64, i64, double, double }* nocapture nonnull readonly align 8 dereferenceable(80) %3, [12 x { i64, i64, {}*, {}*, {}* }]* nocapture nonnull readonly align 8 dereferenceable(480) %4, [14 x { i64, i64, {}*, {}*, {}* }]* nocapture nonnull readonly align 8 dereferenceable(560) %5, { { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, {}*, {}*, {}*, {}*, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* } }* nocapture nonnull readonly align 8 dereferenceable(752) %6) #0 {
pass.3:
  %7 = alloca [3 x {}*], align 8
  %gcframe1625 = alloca [11 x {}*], align 16
  %gcframe1625.sub = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 0
  %.sub = getelementptr inbounds [3 x {}*], [3 x {}*]* %7, i64 0, i64 0
  %8 = bitcast [11 x {}*]* %gcframe1625 to i8*
  call void @llvm.memset.p0i8.i32(i8* noundef nonnull align 16 dereferenceable(88) %8, i8 0, i32 88, i1 false)
  %9 = call {}*** inttoptr (i64 7135162620 to {}*** (i64)*)(i64 260) #6
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:3 within `advance_check2`
; ┌ @ Base.jl:38 within `getproperty`
   %10 = bitcast [11 x {}*]* %gcframe1625 to i64*
   store i64 36, i64* %10, align 16
   %11 = load {}**, {}*** %9, align 8
   %12 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 1
   %13 = bitcast {}** %12 to {}***
   store {}** %11, {}*** %13, align 8
   %14 = bitcast {}*** %9 to {}***
   store {}** %gcframe1625.sub, {}*** %14, align 8
   %15 = getelementptr inbounds { double, double, double, double, double, double, double, double, {}*, {}* }, { double, double, double, double, double, double, double, double, {}*, {}* }* %2, i64 0, i32 0
; └
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:6 within `advance_check2`
; ┌ @ array.jl:126 within `vect`
; │┌ @ array.jl:679 within `_array_for` @ array.jl:676
; ││┌ @ abstractarray.jl:840 within `similar` @ abstractarray.jl:841
; │││┌ @ boot.jl:468 within `Array` @ boot.jl:459
      %16 = call nonnull {}* inttoptr (i64 4320593348 to {}* ({}*, i64)*)({}* inttoptr (i64 4708472368 to {}*), i64 4)
      %17 = bitcast {}* %16 to double**
      %18 = load double*, double** %17, align 8
; │└└└
; │┌ @ array.jl:966 within `setindex!`
    %19 = bitcast double* %18 to <2 x double>*
    store <2 x double> <double 0x3FC5555555555555, double 0x3FD5555555555555>, <2 x double>* %19, align 8
    %20 = getelementptr inbounds double, double* %18, i64 2
    %21 = bitcast double* %20 to <2 x double>*
    store <2 x double> <double 0x3FD5555555555555, double 0x3FC5555555555555>, <2 x double>* %21, align 8
    %22 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 6
    store {}* %16, {}** %22, align 16
; └└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:7 within `advance_check2`
  store {}* inttoptr (i64 4788450624 to {}*), {}** %.sub, align 8
  %23 = getelementptr inbounds [3 x {}*], [3 x {}*]* %7, i64 0, i64 1
  store {}* inttoptr (i64 4788450624 to {}*), {}** %23, align 8
  %24 = getelementptr inbounds [3 x {}*], [3 x {}*]* %7, i64 0, i64 2
  store {}* inttoptr (i64 4317429856 to {}*), {}** %24, align 8
  %25 = call nonnull {}* @j1_vect_2814({}* inttoptr (i64 4721803856 to {}*), {}** nonnull %.sub, i32 3)
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:9 within `advance_check2`
; ┌ @ Base.jl:38 within `getproperty`
   %26 = bitcast {}* %0 to {}**
   %27 = load atomic {}*, {}** %26 unordered, align 8
   %28 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 2
   store {}* %27, {}** %28, align 16
   %29 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 7
   store {}* %25, {}** %29, align 8
; └
; ┌ @ array.jl:369 within `copy`
   %30 = call nonnull {}* inttoptr (i64 4320602180 to {}* ({}*)*)({}* nonnull %27)
; └
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:10 within `advance_check2`
; ┌ @ Base.jl:38 within `getproperty`
   %31 = bitcast {}* %0 to i8*
   %32 = getelementptr inbounds i8, i8* %31, i64 8
   %33 = bitcast i8* %32 to {}**
   %34 = load atomic {}*, {}** %33 unordered, align 8
   %35 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 10
   store {}* %30, {}** %35, align 16
   store {}* %34, {}** %28, align 16
; └
; ┌ @ array.jl:369 within `copy`
   %36 = call nonnull {}* inttoptr (i64 4320602180 to {}* ({}*)*)({}* nonnull %34)
; └
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:11 within `advance_check2`
; ┌ @ Base.jl:38 within `getproperty`
   %37 = getelementptr inbounds i8, i8* %31, i64 16
   %38 = bitcast i8* %37 to {}**
   %39 = load atomic {}*, {}** %38 unordered, align 8
   store {}* %39, {}** %28, align 16
   %40 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 8
   store {}* %36, {}** %40, align 16
; └
; ┌ @ array.jl:369 within `copy`
   %41 = call nonnull {}* inttoptr (i64 4320602180 to {}* ({}*)*)({}* nonnull %39)
; └
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:13 within `advance_check2`
; ┌ @ Base.jl:41 within `dotgetproperty`
; │┌ @ Base.jl:38 within `getproperty`
    %42 = getelementptr inbounds { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }, { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }* %1, i64 0, i32 4
    %43 = load atomic {}*, {}** %42 unordered, align 8
; └└
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ abstractarray.jl:95 within `axes`
; ││┌ @ array.jl:151 within `size`
     %44 = bitcast {}* %43 to { i8*, i64, i16, i16, i32 }*
     %45 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %44, i64 0, i32 1
     %46 = load i64, i64* %45, align 8
; │└└
; │┌ @ broadcast.jl:284 within `instantiate`
; ││┌ @ broadcast.jl:543 within `check_broadcast_axes`
; │││┌ @ abstractarray.jl:95 within `axes`
; ││││┌ @ array.jl:151 within `size`
       %47 = bitcast {}* %30 to { i8*, i64, i16, i16, i32 }*
       %48 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %47, i64 0, i32 1
       %49 = load i64, i64* %48, align 8
; │││└└
; │││┌ @ broadcast.jl:540 within `check_broadcast_shape`
; ││││┌ @ broadcast.jl:518 within `_bcsm`
; │││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
        %50 = icmp ne i64 %46, %49
; │││││└
       %51 = icmp ne i64 %49, 1
; ││││└
      %52 = and i1 %50, %51
      br i1 %52, label %L50, label %L69

L50:                                              ; preds = %pass.3
      %ptls_field1659 = getelementptr inbounds {}**, {}*** %9, i64 2
      %53 = bitcast {}*** %ptls_field1659 to i8**
      %ptls_load16601661 = load i8*, i8** %53, align 8
      %54 = call noalias nonnull {}* @ijl_gc_pool_alloc(i8* %ptls_load16601661, i32 1392, i32 16) #7
      %55 = bitcast {}* %54 to i64*
      %56 = getelementptr inbounds i64, i64* %55, i64 -1
      store atomic i64 4712015424, i64* %56 unordered, align 8
      %57 = bitcast {}* %54 to {}**
      store {}* inttoptr (i64 4725386192 to {}*), {}** %57, align 8
      call void @ijl_throw({}* %54)
      unreachable

L69:                                              ; preds = %pass.3
; │└└└
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:953
; ││┌ @ tuple.jl:398 within `==`
; │││┌ @ tuple.jl:402 within `_eq`
; ││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
       %.not272.not = icmp eq i64 %46, %49
; ││└└└
    br i1 %.not272.not, label %L85, label %L88

L85:                                              ; preds = %L69
    %58 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 9
    store {}* %41, {}** %58, align 8
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954 @ array.jl:346 @ array.jl:322
    %59 = call nonnull {}* @"j__copyto_impl!_2815"({}* nonnull %43, i64 signext 1, {}* nonnull %30, i64 signext 1, i64 signext %46) #0
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954
    br label %L183

L88:                                              ; preds = %L69
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:934 within `broadcast_unalias`
        %.not379 = icmp eq {}* %43, %30
        br i1 %.not379, label %L118, label %L91

L91:                                              ; preds = %L88
; ││││││┌ @ abstractarray.jl:1427 within `unalias`
; │││││││┌ @ abstractarray.jl:1462 within `mightalias`
; ││││││││┌ @ reflection.jl:593 within `isbits`
; │││││││││┌ @ Base.jl:33 within `getproperty`
            %60 = load i8, i8* inttoptr (i64 4708472424 to i8*), align 8
; ││││││││└└
          %61 = and i8 %60, 8
          %.not383.not = icmp eq i8 %61, 0
          br i1 %.not383.not, label %L101, label %L118

L101:                                             ; preds = %L91
; ││││││││┌ @ abstractarray.jl:1486 within `dataids`
; │││││││││┌ @ abstractarray.jl:1187 within `pointer`
; ││││││││││┌ @ pointer.jl:65 within `unsafe_convert`
             %62 = bitcast {}* %43 to i8**
             %63 = load i8*, i8** %62, align 8
             %64 = bitcast {}* %30 to i8**
             %65 = load i8*, i8** %64, align 8
; ││││││││└└└
; ││││││││┌ @ abstractarray.jl:1469 within `_isdisjoint`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %66 = icmp eq i8* %63, %65
; │││││││└└└└
         br i1 %66, label %L113, label %L118

L113:                                             ; preds = %L101
         %67 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 9
         store {}* %41, {}** %67, align 8
; │││││││┌ @ abstractarray.jl:1443 within `unaliascopy`
; ││││││││┌ @ array.jl:369 within `copy`
           %68 = call nonnull {}* inttoptr (i64 4320602180 to {}* ({}*)*)({}* nonnull %30)
; └└└└└└└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:39 within `advance_check2`
; ┌ @ reducedim.jl:1007 within `all`
; │┌ @ reducedim.jl:1007 within `#all#795`
; ││┌ @ reduce.jl:1260 within `_all`
; │││┌ @ array.jl:898 within `iterate`
      br label %L118

L118:                                             ; preds = %L113, %L101, %L91, %L88
      %value_phi198 = phi {}* [ %30, %L88 ], [ %68, %L113 ], [ %30, %L101 ], [ %30, %L91 ]
; └└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:13 within `advance_check2`
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:72 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %.not380 = icmp eq i64 %46, 0
; │││└
     br i1 %.not380, label %L183, label %L139.lr.ph

L139.lr.ph:                                       ; preds = %L118
; ││└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:637 within `extrude`
; ││││││┌ @ broadcast.jl:586 within `newindexer`
; │││││││┌ @ abstractarray.jl:95 within `axes`
; ││││││││┌ @ array.jl:151 within `size`
           %69 = bitcast {}* %value_phi198 to { i8*, i64, i16, i16, i32 }*
           %70 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %69, i64 0, i32 1
           %71 = load i64, i64* %70, align 8
; │││││││└└
; │││││││┌ @ broadcast.jl:587 within `shapeindexer`
; ││││││││┌ @ broadcast.jl:592 within `_newindexer`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %.not382 = icmp eq i64 %71, 1
             %72 = bitcast {}* %value_phi198 to double**
             %73 = load double*, double** %72, align 8
             %74 = bitcast {}* %43 to double**
             %75 = load double*, double** %74, align 8
; ││└└└└└└└└└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:75 within `macro expansion`
     %min.iters.check1086 = icmp ult i64 %46, 4
     br i1 %.not382, label %L139.us.preheader, label %L139.preheader

L139.preheader:                                   ; preds = %L139.lr.ph
     br i1 %min.iters.check1086, label %L139, label %vector.memcheck

vector.memcheck:                                  ; preds = %L139.preheader
     %scevgep = getelementptr double, double* %75, i64 %46
     %scevgep1071 = getelementptr double, double* %73, i64 %46
     %bound0 = icmp ult double* %75, %scevgep1071
     %bound1 = icmp ult double* %73, %scevgep
     %found.conflict = and i1 %bound0, %bound1
     br i1 %found.conflict, label %L139, label %vector.ph

vector.ph:                                        ; preds = %vector.memcheck
     %n.vec = and i64 %46, 9223372036854775804
     br label %vector.body

vector.body:                                      ; preds = %vector.body, %vector.ph
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index = phi i64 [ 0, %vector.ph ], [ %index.next, %vector.body ]
      %76 = getelementptr inbounds double, double* %73, i64 %index
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %77 = bitcast double* %76 to <2 x double>*
          %wide.load = load <2 x double>, <2 x double>* %77, align 8
          %78 = getelementptr inbounds double, double* %76, i64 2
          %79 = bitcast double* %78 to <2 x double>*
          %wide.load1073 = load <2 x double>, <2 x double>* %79, align 8
; │││└└└└└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %80 = getelementptr inbounds double, double* %75, i64 %index
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ array.jl:966 within `setindex!`
      %81 = bitcast double* %80 to <2 x double>*
      store <2 x double> %wide.load, <2 x double>* %81, align 8
      %82 = getelementptr inbounds double, double* %80, i64 2
      %83 = bitcast double* %82 to <2 x double>*
      store <2 x double> %wide.load1073, <2 x double>* %83, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index.next = add nuw i64 %index, 4
      %84 = icmp eq i64 %index.next, %n.vec
      br i1 %84, label %middle.block, label %vector.body

middle.block:                                     ; preds = %vector.body
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n = icmp eq i64 %46, %n.vec
     br i1 %cmp.n, label %L183, label %L139

L139.us.preheader:                                ; preds = %L139.lr.ph
     br i1 %min.iters.check1086, label %L139.us, label %vector.memcheck1074

vector.memcheck1074:                              ; preds = %L139.us.preheader
     %scevgep1075 = getelementptr double, double* %75, i64 %46
     %scevgep1077 = getelementptr double, double* %73, i64 1
     %bound01079 = icmp ult double* %75, %scevgep1077
     %bound11080 = icmp ult double* %73, %scevgep1075
     %found.conflict1081 = and i1 %bound01079, %bound11080
     br i1 %found.conflict1081, label %L139.us, label %vector.ph1087

vector.ph1087:                                    ; preds = %vector.memcheck1074
     %n.vec1089 = and i64 %46, 9223372036854775804
     br label %vector.body1085

vector.body1085:                                  ; preds = %vector.body1085, %vector.ph1087
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1090 = phi i64 [ 0, %vector.ph1087 ], [ %index.next1091, %vector.body1085 ]
      %85 = load double, double* %73, align 8
      %broadcast.splatinsert = insertelement <2 x double> poison, double %85, i32 0
      %broadcast.splat = shufflevector <2 x double> %broadcast.splatinsert, <2 x double> poison, <2 x i32> zeroinitializer
      %86 = getelementptr inbounds double, double* %75, i64 %index1090
      %87 = bitcast double* %86 to <2 x double>*
      store <2 x double> %broadcast.splat, <2 x double>* %87, align 8
      %88 = getelementptr inbounds double, double* %86, i64 2
      %89 = bitcast double* %88 to <2 x double>*
      store <2 x double> %broadcast.splat, <2 x double>* %89, align 8
      %index.next1091 = add nuw i64 %index1090, 4
      %90 = icmp eq i64 %index.next1091, %n.vec1089
      br i1 %90, label %middle.block1083, label %vector.body1085

middle.block1083:                                 ; preds = %vector.body1085
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1093 = icmp eq i64 %46, %n.vec1089
     br i1 %cmp.n1093, label %L183, label %L139.us

L139.us:                                          ; preds = %L139.us, %middle.block1083, %vector.memcheck1074, %L139.us.preheader
     %value_phi199479.us = phi i64 [ %93, %L139.us ], [ %n.vec1089, %middle.block1083 ], [ 0, %L139.us.preheader ], [ 0, %vector.memcheck1074 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %91 = load double, double* %73, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %92 = getelementptr inbounds double, double* %75, i64 %value_phi199479.us
      store double %91, double* %92, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %93 = add nuw nsw i64 %value_phi199479.us, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond849.not = icmp eq i64 %93, %46
; │││└
     br i1 %exitcond849.not, label %L183, label %L139.us

L139:                                             ; preds = %L139, %middle.block, %vector.memcheck, %L139.preheader
     %value_phi199479 = phi i64 [ %97, %L139 ], [ %n.vec, %middle.block ], [ 0, %L139.preheader ], [ 0, %vector.memcheck ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %94 = getelementptr inbounds double, double* %73, i64 %value_phi199479
          %95 = load double, double* %94, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %96 = getelementptr inbounds double, double* %75, i64 %value_phi199479
      store double %95, double* %96, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %97 = add nuw nsw i64 %value_phi199479, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond843.not = icmp eq i64 %97, %46
; │││└
     br i1 %exitcond843.not, label %L183, label %L139

L183:                                             ; preds = %L139, %L139.us, %middle.block1083, %middle.block, %L118, %L85
; └└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:14 within `advance_check2`
; ┌ @ Base.jl:41 within `dotgetproperty`
; │┌ @ Base.jl:38 within `getproperty`
    %98 = getelementptr inbounds { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }, { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }* %1, i64 0, i32 5
    %99 = load atomic {}*, {}** %98 unordered, align 8
; └└
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ abstractarray.jl:95 within `axes`
; ││┌ @ array.jl:151 within `size`
     %100 = bitcast {}* %99 to { i8*, i64, i16, i16, i32 }*
     %101 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %100, i64 0, i32 1
     %102 = load i64, i64* %101, align 8
; │└└
; │┌ @ broadcast.jl:284 within `instantiate`
; ││┌ @ broadcast.jl:543 within `check_broadcast_axes`
; │││┌ @ abstractarray.jl:95 within `axes`
; ││││┌ @ array.jl:151 within `size`
       %103 = bitcast {}* %36 to { i8*, i64, i16, i16, i32 }*
       %104 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %103, i64 0, i32 1
       %105 = load i64, i64* %104, align 8
; │││└└
; │││┌ @ broadcast.jl:540 within `check_broadcast_shape`
; ││││┌ @ broadcast.jl:518 within `_bcsm`
; │││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
        %106 = icmp ne i64 %102, %105
; │││││└
       %107 = icmp ne i64 %105, 1
; ││││└
      %108 = and i1 %106, %107
      br i1 %108, label %L202, label %L221

L202:                                             ; preds = %L183
      %ptls_field16031656 = getelementptr inbounds {}**, {}*** %9, i64 2
      %109 = bitcast {}*** %ptls_field16031656 to i8**
      %ptls_load160416571658 = load i8*, i8** %109, align 8
      %110 = call noalias nonnull {}* @ijl_gc_pool_alloc(i8* %ptls_load160416571658, i32 1392, i32 16) #7
      %111 = bitcast {}* %110 to i64*
      %112 = getelementptr inbounds i64, i64* %111, i64 -1
      store atomic i64 4712015424, i64* %112 unordered, align 8
      %113 = bitcast {}* %110 to {}**
      store {}* inttoptr (i64 4725386192 to {}*), {}** %113, align 8
      call void @ijl_throw({}* %110)
      unreachable

L221:                                             ; preds = %L183
; │└└└
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:953
; ││┌ @ tuple.jl:398 within `==`
; │││┌ @ tuple.jl:402 within `_eq`
; ││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
       %.not274.not = icmp eq i64 %102, %105
; ││└└└
    br i1 %.not274.not, label %L237, label %L240

L237:                                             ; preds = %L221
    %114 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 9
    store {}* %41, {}** %114, align 8
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954 @ array.jl:346 @ array.jl:322
    %115 = call nonnull {}* @"j__copyto_impl!_2816"({}* nonnull %99, i64 signext 1, {}* nonnull %36, i64 signext 1, i64 signext %102) #0
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954
    br label %L335

L240:                                             ; preds = %L221
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:934 within `broadcast_unalias`
        %.not372 = icmp eq {}* %99, %36
        br i1 %.not372, label %L270, label %L243

L243:                                             ; preds = %L240
; ││││││┌ @ abstractarray.jl:1427 within `unalias`
; │││││││┌ @ abstractarray.jl:1462 within `mightalias`
; ││││││││┌ @ reflection.jl:593 within `isbits`
; │││││││││┌ @ Base.jl:33 within `getproperty`
            %116 = load i8, i8* inttoptr (i64 4708472424 to i8*), align 8
; ││││││││└└
          %117 = and i8 %116, 8
          %.not376.not = icmp eq i8 %117, 0
          br i1 %.not376.not, label %L253, label %L270

L253:                                             ; preds = %L243
; ││││││││┌ @ abstractarray.jl:1486 within `dataids`
; │││││││││┌ @ abstractarray.jl:1187 within `pointer`
; ││││││││││┌ @ pointer.jl:65 within `unsafe_convert`
             %118 = bitcast {}* %99 to i8**
             %119 = load i8*, i8** %118, align 8
             %120 = bitcast {}* %36 to i8**
             %121 = load i8*, i8** %120, align 8
; ││││││││└└└
; ││││││││┌ @ abstractarray.jl:1469 within `_isdisjoint`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %122 = icmp eq i8* %119, %121
; │││││││└└└└
         br i1 %122, label %L265, label %L270

L265:                                             ; preds = %L253
         %123 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 9
         store {}* %41, {}** %123, align 8
; │││││││┌ @ abstractarray.jl:1443 within `unaliascopy`
; ││││││││┌ @ array.jl:369 within `copy`
           %124 = call nonnull {}* inttoptr (i64 4320602180 to {}* ({}*)*)({}* nonnull %36)
; └└└└└└└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:39 within `advance_check2`
; ┌ @ reducedim.jl:1007 within `all`
; │┌ @ reducedim.jl:1007 within `#all#795`
; ││┌ @ reduce.jl:1260 within `_all`
; │││┌ @ array.jl:898 within `iterate`
      br label %L270

L270:                                             ; preds = %L265, %L253, %L243, %L240
      %value_phi188 = phi {}* [ %36, %L240 ], [ %124, %L265 ], [ %36, %L253 ], [ %36, %L243 ]
; └└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:14 within `advance_check2`
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:72 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %.not373 = icmp eq i64 %102, 0
; │││└
     br i1 %.not373, label %L335, label %L291.lr.ph

L291.lr.ph:                                       ; preds = %L270
; ││└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:637 within `extrude`
; ││││││┌ @ broadcast.jl:586 within `newindexer`
; │││││││┌ @ abstractarray.jl:95 within `axes`
; ││││││││┌ @ array.jl:151 within `size`
           %125 = bitcast {}* %value_phi188 to { i8*, i64, i16, i16, i32 }*
           %126 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %125, i64 0, i32 1
           %127 = load i64, i64* %126, align 8
; │││││││└└
; │││││││┌ @ broadcast.jl:587 within `shapeindexer`
; ││││││││┌ @ broadcast.jl:592 within `_newindexer`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %.not375 = icmp eq i64 %127, 1
             %128 = bitcast {}* %value_phi188 to double**
             %129 = load double*, double** %128, align 8
             %130 = bitcast {}* %99 to double**
             %131 = load double*, double** %130, align 8
; ││└└└└└└└└└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:75 within `macro expansion`
     %min.iters.check1130 = icmp ult i64 %102, 4
     br i1 %.not375, label %L291.us.preheader, label %L291.preheader

L291.preheader:                                   ; preds = %L291.lr.ph
     br i1 %min.iters.check1130, label %L291, label %vector.memcheck1096

vector.memcheck1096:                              ; preds = %L291.preheader
     %scevgep1097 = getelementptr double, double* %131, i64 %102
     %scevgep1099 = getelementptr double, double* %129, i64 %102
     %bound01101 = icmp ult double* %131, %scevgep1099
     %bound11102 = icmp ult double* %129, %scevgep1097
     %found.conflict1103 = and i1 %bound01101, %bound11102
     br i1 %found.conflict1103, label %L291, label %vector.ph1109

vector.ph1109:                                    ; preds = %vector.memcheck1096
     %n.vec1111 = and i64 %102, 9223372036854775804
     br label %vector.body1107

vector.body1107:                                  ; preds = %vector.body1107, %vector.ph1109
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1112 = phi i64 [ 0, %vector.ph1109 ], [ %index.next1113, %vector.body1107 ]
      %132 = getelementptr inbounds double, double* %129, i64 %index1112
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %133 = bitcast double* %132 to <2 x double>*
          %wide.load1116 = load <2 x double>, <2 x double>* %133, align 8
          %134 = getelementptr inbounds double, double* %132, i64 2
          %135 = bitcast double* %134 to <2 x double>*
          %wide.load1117 = load <2 x double>, <2 x double>* %135, align 8
; │││└└└└└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %136 = getelementptr inbounds double, double* %131, i64 %index1112
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ array.jl:966 within `setindex!`
      %137 = bitcast double* %136 to <2 x double>*
      store <2 x double> %wide.load1116, <2 x double>* %137, align 8
      %138 = getelementptr inbounds double, double* %136, i64 2
      %139 = bitcast double* %138 to <2 x double>*
      store <2 x double> %wide.load1117, <2 x double>* %139, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index.next1113 = add nuw i64 %index1112, 4
      %140 = icmp eq i64 %index.next1113, %n.vec1111
      br i1 %140, label %middle.block1105, label %vector.body1107

middle.block1105:                                 ; preds = %vector.body1107
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1115 = icmp eq i64 %102, %n.vec1111
     br i1 %cmp.n1115, label %L335, label %L291

L291.us.preheader:                                ; preds = %L291.lr.ph
     br i1 %min.iters.check1130, label %L291.us, label %vector.memcheck1118

vector.memcheck1118:                              ; preds = %L291.us.preheader
     %scevgep1119 = getelementptr double, double* %131, i64 %102
     %scevgep1121 = getelementptr double, double* %129, i64 1
     %bound01123 = icmp ult double* %131, %scevgep1121
     %bound11124 = icmp ult double* %129, %scevgep1119
     %found.conflict1125 = and i1 %bound01123, %bound11124
     br i1 %found.conflict1125, label %L291.us, label %vector.ph1131

vector.ph1131:                                    ; preds = %vector.memcheck1118
     %n.vec1133 = and i64 %102, 9223372036854775804
     br label %vector.body1129

vector.body1129:                                  ; preds = %vector.body1129, %vector.ph1131
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1134 = phi i64 [ 0, %vector.ph1131 ], [ %index.next1135, %vector.body1129 ]
      %141 = load double, double* %129, align 8
      %broadcast.splatinsert1138 = insertelement <2 x double> poison, double %141, i32 0
      %broadcast.splat1139 = shufflevector <2 x double> %broadcast.splatinsert1138, <2 x double> poison, <2 x i32> zeroinitializer
      %142 = getelementptr inbounds double, double* %131, i64 %index1134
      %143 = bitcast double* %142 to <2 x double>*
      store <2 x double> %broadcast.splat1139, <2 x double>* %143, align 8
      %144 = getelementptr inbounds double, double* %142, i64 2
      %145 = bitcast double* %144 to <2 x double>*
      store <2 x double> %broadcast.splat1139, <2 x double>* %145, align 8
      %index.next1135 = add nuw i64 %index1134, 4
      %146 = icmp eq i64 %index.next1135, %n.vec1133
      br i1 %146, label %middle.block1127, label %vector.body1129

middle.block1127:                                 ; preds = %vector.body1129
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1137 = icmp eq i64 %102, %n.vec1133
     br i1 %cmp.n1137, label %L335, label %L291.us

L291.us:                                          ; preds = %L291.us, %middle.block1127, %vector.memcheck1118, %L291.us.preheader
     %value_phi189477.us = phi i64 [ %149, %L291.us ], [ %n.vec1133, %middle.block1127 ], [ 0, %L291.us.preheader ], [ 0, %vector.memcheck1118 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %147 = load double, double* %129, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %148 = getelementptr inbounds double, double* %131, i64 %value_phi189477.us
      store double %147, double* %148, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %149 = add nuw nsw i64 %value_phi189477.us, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond848.not = icmp eq i64 %149, %102
; │││└
     br i1 %exitcond848.not, label %L335, label %L291.us

L291:                                             ; preds = %L291, %middle.block1105, %vector.memcheck1096, %L291.preheader
     %value_phi189477 = phi i64 [ %153, %L291 ], [ %n.vec1111, %middle.block1105 ], [ 0, %L291.preheader ], [ 0, %vector.memcheck1096 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %150 = getelementptr inbounds double, double* %129, i64 %value_phi189477
          %151 = load double, double* %150, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %152 = getelementptr inbounds double, double* %131, i64 %value_phi189477
      store double %151, double* %152, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %153 = add nuw nsw i64 %value_phi189477, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond842.not = icmp eq i64 %153, %102
; │││└
     br i1 %exitcond842.not, label %L335, label %L291

L335:                                             ; preds = %L291, %L291.us, %middle.block1127, %middle.block1105, %L270, %L237
; └└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:15 within `advance_check2`
; ┌ @ Base.jl:41 within `dotgetproperty`
; │┌ @ Base.jl:38 within `getproperty`
    %154 = getelementptr inbounds { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }, { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }* %1, i64 0, i32 6
    %155 = load atomic {}*, {}** %154 unordered, align 8
; └└
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ abstractarray.jl:95 within `axes`
; ││┌ @ array.jl:151 within `size`
     %156 = bitcast {}* %155 to { i8*, i64, i16, i16, i32 }*
     %157 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %156, i64 0, i32 1
     %158 = load i64, i64* %157, align 8
; │└└
; │┌ @ broadcast.jl:284 within `instantiate`
; ││┌ @ broadcast.jl:543 within `check_broadcast_axes`
; │││┌ @ abstractarray.jl:95 within `axes`
; ││││┌ @ array.jl:151 within `size`
       %159 = bitcast {}* %41 to { i8*, i64, i16, i16, i32 }*
       %160 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %159, i64 0, i32 1
       %161 = load i64, i64* %160, align 8
; │││└└
; │││┌ @ broadcast.jl:540 within `check_broadcast_shape`
; ││││┌ @ broadcast.jl:518 within `_bcsm`
; │││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
        %162 = icmp ne i64 %158, %161
; │││││└
       %163 = icmp ne i64 %161, 1
; ││││└
      %164 = and i1 %162, %163
      br i1 %164, label %L354, label %L373

L354:                                             ; preds = %L335
      %ptls_field16051653 = getelementptr inbounds {}**, {}*** %9, i64 2
      %165 = bitcast {}*** %ptls_field16051653 to i8**
      %ptls_load160616541655 = load i8*, i8** %165, align 8
      %166 = call noalias nonnull {}* @ijl_gc_pool_alloc(i8* %ptls_load160616541655, i32 1392, i32 16) #7
      %167 = bitcast {}* %166 to i64*
      %168 = getelementptr inbounds i64, i64* %167, i64 -1
      store atomic i64 4712015424, i64* %168 unordered, align 8
      %169 = bitcast {}* %166 to {}**
      store {}* inttoptr (i64 4725386192 to {}*), {}** %169, align 8
      call void @ijl_throw({}* %166)
      unreachable

L373:                                             ; preds = %L335
; │└└└
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:953
; ││┌ @ tuple.jl:398 within `==`
; │││┌ @ tuple.jl:402 within `_eq`
; ││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
       %.not276.not = icmp eq i64 %158, %161
; ││└└└
    br i1 %.not276.not, label %L389, label %L392

L389:                                             ; preds = %L373
    %170 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 9
    store {}* %41, {}** %170, align 8
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954 @ array.jl:346 @ array.jl:322
    %171 = call nonnull {}* @"j__copyto_impl!_2817"({}* nonnull %155, i64 signext 1, {}* nonnull %41, i64 signext 1, i64 signext %158) #0
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954
    br label %L487

L392:                                             ; preds = %L373
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:934 within `broadcast_unalias`
        %.not365 = icmp eq {}* %155, %41
        br i1 %.not365, label %L422, label %L395

L395:                                             ; preds = %L392
; ││││││┌ @ abstractarray.jl:1427 within `unalias`
; │││││││┌ @ abstractarray.jl:1462 within `mightalias`
; ││││││││┌ @ reflection.jl:593 within `isbits`
; │││││││││┌ @ Base.jl:33 within `getproperty`
            %172 = load i8, i8* inttoptr (i64 4708472424 to i8*), align 8
; ││││││││└└
          %173 = and i8 %172, 8
          %.not369.not = icmp eq i8 %173, 0
          br i1 %.not369.not, label %L405, label %L422

L405:                                             ; preds = %L395
; ││││││││┌ @ abstractarray.jl:1486 within `dataids`
; │││││││││┌ @ abstractarray.jl:1187 within `pointer`
; ││││││││││┌ @ pointer.jl:65 within `unsafe_convert`
             %174 = bitcast {}* %155 to i8**
             %175 = load i8*, i8** %174, align 8
             %176 = bitcast {}* %41 to i8**
             %177 = load i8*, i8** %176, align 8
; ││││││││└└└
; ││││││││┌ @ abstractarray.jl:1469 within `_isdisjoint`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %178 = icmp eq i8* %175, %177
; │││││││└└└└
         br i1 %178, label %L417, label %L422

L417:                                             ; preds = %L405
         %179 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 9
         store {}* %41, {}** %179, align 8
; │││││││┌ @ abstractarray.jl:1443 within `unaliascopy`
; ││││││││┌ @ array.jl:369 within `copy`
           %180 = call nonnull {}* inttoptr (i64 4320602180 to {}* ({}*)*)({}* nonnull %41)
; └└└└└└└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:39 within `advance_check2`
; ┌ @ reducedim.jl:1007 within `all`
; │┌ @ reducedim.jl:1007 within `#all#795`
; ││┌ @ reduce.jl:1260 within `_all`
; │││┌ @ array.jl:898 within `iterate`
      br label %L422

L422:                                             ; preds = %L417, %L405, %L395, %L392
      %value_phi178 = phi {}* [ %41, %L392 ], [ %180, %L417 ], [ %41, %L405 ], [ %41, %L395 ]
; └└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:15 within `advance_check2`
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:72 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %.not366 = icmp eq i64 %158, 0
; │││└
     br i1 %.not366, label %L487, label %L443.lr.ph

L443.lr.ph:                                       ; preds = %L422
; ││└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:637 within `extrude`
; ││││││┌ @ broadcast.jl:586 within `newindexer`
; │││││││┌ @ abstractarray.jl:95 within `axes`
; ││││││││┌ @ array.jl:151 within `size`
           %181 = bitcast {}* %value_phi178 to { i8*, i64, i16, i16, i32 }*
           %182 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %181, i64 0, i32 1
           %183 = load i64, i64* %182, align 8
; │││││││└└
; │││││││┌ @ broadcast.jl:587 within `shapeindexer`
; ││││││││┌ @ broadcast.jl:592 within `_newindexer`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %.not368 = icmp eq i64 %183, 1
             %184 = bitcast {}* %value_phi178 to double**
             %185 = load double*, double** %184, align 8
             %186 = bitcast {}* %155 to double**
             %187 = load double*, double** %186, align 8
; ││└└└└└└└└└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:75 within `macro expansion`
     %min.iters.check1176 = icmp ult i64 %158, 4
     br i1 %.not368, label %L443.us.preheader, label %L443.preheader

L443.preheader:                                   ; preds = %L443.lr.ph
     br i1 %min.iters.check1176, label %L443, label %vector.memcheck1142

vector.memcheck1142:                              ; preds = %L443.preheader
     %scevgep1143 = getelementptr double, double* %187, i64 %158
     %scevgep1145 = getelementptr double, double* %185, i64 %158
     %bound01147 = icmp ult double* %187, %scevgep1145
     %bound11148 = icmp ult double* %185, %scevgep1143
     %found.conflict1149 = and i1 %bound01147, %bound11148
     br i1 %found.conflict1149, label %L443, label %vector.ph1155

vector.ph1155:                                    ; preds = %vector.memcheck1142
     %n.vec1157 = and i64 %158, 9223372036854775804
     br label %vector.body1153

vector.body1153:                                  ; preds = %vector.body1153, %vector.ph1155
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1158 = phi i64 [ 0, %vector.ph1155 ], [ %index.next1159, %vector.body1153 ]
      %188 = getelementptr inbounds double, double* %185, i64 %index1158
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %189 = bitcast double* %188 to <2 x double>*
          %wide.load1162 = load <2 x double>, <2 x double>* %189, align 8
          %190 = getelementptr inbounds double, double* %188, i64 2
          %191 = bitcast double* %190 to <2 x double>*
          %wide.load1163 = load <2 x double>, <2 x double>* %191, align 8
; │││└└└└└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %192 = getelementptr inbounds double, double* %187, i64 %index1158
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ array.jl:966 within `setindex!`
      %193 = bitcast double* %192 to <2 x double>*
      store <2 x double> %wide.load1162, <2 x double>* %193, align 8
      %194 = getelementptr inbounds double, double* %192, i64 2
      %195 = bitcast double* %194 to <2 x double>*
      store <2 x double> %wide.load1163, <2 x double>* %195, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index.next1159 = add nuw i64 %index1158, 4
      %196 = icmp eq i64 %index.next1159, %n.vec1157
      br i1 %196, label %middle.block1151, label %vector.body1153

middle.block1151:                                 ; preds = %vector.body1153
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1161 = icmp eq i64 %158, %n.vec1157
     br i1 %cmp.n1161, label %L487, label %L443

L443.us.preheader:                                ; preds = %L443.lr.ph
     br i1 %min.iters.check1176, label %L443.us, label %vector.memcheck1164

vector.memcheck1164:                              ; preds = %L443.us.preheader
     %scevgep1165 = getelementptr double, double* %187, i64 %158
     %scevgep1167 = getelementptr double, double* %185, i64 1
     %bound01169 = icmp ult double* %187, %scevgep1167
     %bound11170 = icmp ult double* %185, %scevgep1165
     %found.conflict1171 = and i1 %bound01169, %bound11170
     br i1 %found.conflict1171, label %L443.us, label %vector.ph1177

vector.ph1177:                                    ; preds = %vector.memcheck1164
     %n.vec1179 = and i64 %158, 9223372036854775804
     br label %vector.body1175

vector.body1175:                                  ; preds = %vector.body1175, %vector.ph1177
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1180 = phi i64 [ 0, %vector.ph1177 ], [ %index.next1181, %vector.body1175 ]
      %197 = load double, double* %185, align 8
      %broadcast.splatinsert1184 = insertelement <2 x double> poison, double %197, i32 0
      %broadcast.splat1185 = shufflevector <2 x double> %broadcast.splatinsert1184, <2 x double> poison, <2 x i32> zeroinitializer
      %198 = getelementptr inbounds double, double* %187, i64 %index1180
      %199 = bitcast double* %198 to <2 x double>*
      store <2 x double> %broadcast.splat1185, <2 x double>* %199, align 8
      %200 = getelementptr inbounds double, double* %198, i64 2
      %201 = bitcast double* %200 to <2 x double>*
      store <2 x double> %broadcast.splat1185, <2 x double>* %201, align 8
      %index.next1181 = add nuw i64 %index1180, 4
      %202 = icmp eq i64 %index.next1181, %n.vec1179
      br i1 %202, label %middle.block1173, label %vector.body1175

middle.block1173:                                 ; preds = %vector.body1175
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1183 = icmp eq i64 %158, %n.vec1179
     br i1 %cmp.n1183, label %L487, label %L443.us

L443.us:                                          ; preds = %L443.us, %middle.block1173, %vector.memcheck1164, %L443.us.preheader
     %value_phi179475.us = phi i64 [ %205, %L443.us ], [ %n.vec1179, %middle.block1173 ], [ 0, %L443.us.preheader ], [ 0, %vector.memcheck1164 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %203 = load double, double* %185, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %204 = getelementptr inbounds double, double* %187, i64 %value_phi179475.us
      store double %203, double* %204, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %205 = add nuw nsw i64 %value_phi179475.us, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond847.not = icmp eq i64 %205, %158
; │││└
     br i1 %exitcond847.not, label %L487, label %L443.us

L443:                                             ; preds = %L443, %middle.block1151, %vector.memcheck1142, %L443.preheader
     %value_phi179475 = phi i64 [ %209, %L443 ], [ %n.vec1157, %middle.block1151 ], [ 0, %L443.preheader ], [ 0, %vector.memcheck1142 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %206 = getelementptr inbounds double, double* %185, i64 %value_phi179475
          %207 = load double, double* %206, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %208 = getelementptr inbounds double, double* %187, i64 %value_phi179475
      store double %207, double* %208, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %209 = add nuw nsw i64 %value_phi179475, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond841.not = icmp eq i64 %209, %158
; │││└
     br i1 %exitcond841.not, label %L487, label %L443

L487:                                             ; preds = %L443, %L443.us, %middle.block1173, %middle.block1151, %L422, %L389
; └└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:17 within `advance_check2`
; ┌ @ Base.jl:41 within `dotgetproperty`
; │┌ @ Base.jl:38 within `getproperty`
    %210 = getelementptr inbounds { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }, { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }* %1, i64 0, i32 7
    %211 = load atomic {}*, {}** %210 unordered, align 8
; └└
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ abstractarray.jl:95 within `axes`
; ││┌ @ array.jl:151 within `size`
     %212 = bitcast {}* %211 to { i8*, i64, i16, i16, i32 }*
     %213 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %212, i64 0, i32 1
     %214 = load i64, i64* %213, align 8
; │└└
; │┌ @ broadcast.jl:284 within `instantiate`
; ││┌ @ broadcast.jl:543 within `check_broadcast_axes`
; │││┌ @ abstractarray.jl:95 within `axes`
; ││││┌ @ array.jl:151 within `size`
       %215 = load i64, i64* %48, align 8
; │││└└
; │││┌ @ broadcast.jl:540 within `check_broadcast_shape`
; ││││┌ @ broadcast.jl:518 within `_bcsm`
; │││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
        %216 = icmp ne i64 %214, %215
; │││││└
       %217 = icmp ne i64 %215, 1
; ││││└
      %218 = and i1 %216, %217
      br i1 %218, label %L506, label %L525

L506:                                             ; preds = %L487
      %ptls_field16071650 = getelementptr inbounds {}**, {}*** %9, i64 2
      %219 = bitcast {}*** %ptls_field16071650 to i8**
      %ptls_load160816511652 = load i8*, i8** %219, align 8
      %220 = call noalias nonnull {}* @ijl_gc_pool_alloc(i8* %ptls_load160816511652, i32 1392, i32 16) #7
      %221 = bitcast {}* %220 to i64*
      %222 = getelementptr inbounds i64, i64* %221, i64 -1
      store atomic i64 4712015424, i64* %222 unordered, align 8
      %223 = bitcast {}* %220 to {}**
      store {}* inttoptr (i64 4725386192 to {}*), {}** %223, align 8
      call void @ijl_throw({}* %220)
      unreachable

L525:                                             ; preds = %L487
; │└└└
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:953
; ││┌ @ tuple.jl:398 within `==`
; │││┌ @ tuple.jl:402 within `_eq`
; ││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
       %.not278.not = icmp eq i64 %214, %215
; ││└└└
    br i1 %.not278.not, label %L541, label %L544

L541:                                             ; preds = %L525
    %224 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 9
    store {}* %41, {}** %224, align 8
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954 @ array.jl:346 @ array.jl:322
    %225 = call nonnull {}* @"j__copyto_impl!_2818"({}* nonnull %211, i64 signext 1, {}* nonnull %30, i64 signext 1, i64 signext %214) #0
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954
    br label %L639

L544:                                             ; preds = %L525
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:934 within `broadcast_unalias`
        %.not358 = icmp eq {}* %211, %30
        br i1 %.not358, label %L574, label %L547

L547:                                             ; preds = %L544
; ││││││┌ @ abstractarray.jl:1427 within `unalias`
; │││││││┌ @ abstractarray.jl:1462 within `mightalias`
; ││││││││┌ @ reflection.jl:593 within `isbits`
; │││││││││┌ @ Base.jl:33 within `getproperty`
            %226 = load i8, i8* inttoptr (i64 4708472424 to i8*), align 8
; ││││││││└└
          %227 = and i8 %226, 8
          %.not362.not = icmp eq i8 %227, 0
          br i1 %.not362.not, label %L557, label %L574

L557:                                             ; preds = %L547
; ││││││││┌ @ abstractarray.jl:1486 within `dataids`
; │││││││││┌ @ abstractarray.jl:1187 within `pointer`
; ││││││││││┌ @ pointer.jl:65 within `unsafe_convert`
             %228 = bitcast {}* %211 to i8**
             %229 = load i8*, i8** %228, align 8
             %230 = bitcast {}* %30 to i8**
             %231 = load i8*, i8** %230, align 8
; ││││││││└└└
; ││││││││┌ @ abstractarray.jl:1469 within `_isdisjoint`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %232 = icmp eq i8* %229, %231
; │││││││└└└└
         br i1 %232, label %L569, label %L574

L569:                                             ; preds = %L557
         %233 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 9
         store {}* %41, {}** %233, align 8
; │││││││┌ @ abstractarray.jl:1443 within `unaliascopy`
; ││││││││┌ @ array.jl:369 within `copy`
           %234 = call nonnull {}* inttoptr (i64 4320602180 to {}* ({}*)*)({}* nonnull %30)
; └└└└└└└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:39 within `advance_check2`
; ┌ @ reducedim.jl:1007 within `all`
; │┌ @ reducedim.jl:1007 within `#all#795`
; ││┌ @ reduce.jl:1260 within `_all`
; │││┌ @ array.jl:898 within `iterate`
      br label %L574

L574:                                             ; preds = %L569, %L557, %L547, %L544
      %value_phi168 = phi {}* [ %30, %L544 ], [ %234, %L569 ], [ %30, %L557 ], [ %30, %L547 ]
; └└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:17 within `advance_check2`
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:72 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %.not359 = icmp eq i64 %214, 0
; │││└
     br i1 %.not359, label %L639, label %L595.lr.ph

L595.lr.ph:                                       ; preds = %L574
; ││└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:637 within `extrude`
; ││││││┌ @ broadcast.jl:586 within `newindexer`
; │││││││┌ @ abstractarray.jl:95 within `axes`
; ││││││││┌ @ array.jl:151 within `size`
           %235 = bitcast {}* %value_phi168 to { i8*, i64, i16, i16, i32 }*
           %236 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %235, i64 0, i32 1
           %237 = load i64, i64* %236, align 8
; │││││││└└
; │││││││┌ @ broadcast.jl:587 within `shapeindexer`
; ││││││││┌ @ broadcast.jl:592 within `_newindexer`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %.not361 = icmp eq i64 %237, 1
             %238 = bitcast {}* %value_phi168 to double**
             %239 = load double*, double** %238, align 8
             %240 = bitcast {}* %211 to double**
             %241 = load double*, double** %240, align 8
; ││└└└└└└└└└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:75 within `macro expansion`
     %min.iters.check1222 = icmp ult i64 %214, 4
     br i1 %.not361, label %L595.us.preheader, label %L595.preheader

L595.preheader:                                   ; preds = %L595.lr.ph
     br i1 %min.iters.check1222, label %L595, label %vector.memcheck1188

vector.memcheck1188:                              ; preds = %L595.preheader
     %scevgep1189 = getelementptr double, double* %241, i64 %214
     %scevgep1191 = getelementptr double, double* %239, i64 %214
     %bound01193 = icmp ult double* %241, %scevgep1191
     %bound11194 = icmp ult double* %239, %scevgep1189
     %found.conflict1195 = and i1 %bound01193, %bound11194
     br i1 %found.conflict1195, label %L595, label %vector.ph1201

vector.ph1201:                                    ; preds = %vector.memcheck1188
     %n.vec1203 = and i64 %214, 9223372036854775804
     br label %vector.body1199

vector.body1199:                                  ; preds = %vector.body1199, %vector.ph1201
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1204 = phi i64 [ 0, %vector.ph1201 ], [ %index.next1205, %vector.body1199 ]
      %242 = getelementptr inbounds double, double* %239, i64 %index1204
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %243 = bitcast double* %242 to <2 x double>*
          %wide.load1208 = load <2 x double>, <2 x double>* %243, align 8
          %244 = getelementptr inbounds double, double* %242, i64 2
          %245 = bitcast double* %244 to <2 x double>*
          %wide.load1209 = load <2 x double>, <2 x double>* %245, align 8
; │││└└└└└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %246 = getelementptr inbounds double, double* %241, i64 %index1204
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ array.jl:966 within `setindex!`
      %247 = bitcast double* %246 to <2 x double>*
      store <2 x double> %wide.load1208, <2 x double>* %247, align 8
      %248 = getelementptr inbounds double, double* %246, i64 2
      %249 = bitcast double* %248 to <2 x double>*
      store <2 x double> %wide.load1209, <2 x double>* %249, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index.next1205 = add nuw i64 %index1204, 4
      %250 = icmp eq i64 %index.next1205, %n.vec1203
      br i1 %250, label %middle.block1197, label %vector.body1199

middle.block1197:                                 ; preds = %vector.body1199
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1207 = icmp eq i64 %214, %n.vec1203
     br i1 %cmp.n1207, label %L639, label %L595

L595.us.preheader:                                ; preds = %L595.lr.ph
     br i1 %min.iters.check1222, label %L595.us, label %vector.memcheck1210

vector.memcheck1210:                              ; preds = %L595.us.preheader
     %scevgep1211 = getelementptr double, double* %241, i64 %214
     %scevgep1213 = getelementptr double, double* %239, i64 1
     %bound01215 = icmp ult double* %241, %scevgep1213
     %bound11216 = icmp ult double* %239, %scevgep1211
     %found.conflict1217 = and i1 %bound01215, %bound11216
     br i1 %found.conflict1217, label %L595.us, label %vector.ph1223

vector.ph1223:                                    ; preds = %vector.memcheck1210
     %n.vec1225 = and i64 %214, 9223372036854775804
     br label %vector.body1221

vector.body1221:                                  ; preds = %vector.body1221, %vector.ph1223
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1226 = phi i64 [ 0, %vector.ph1223 ], [ %index.next1227, %vector.body1221 ]
      %251 = load double, double* %239, align 8
      %broadcast.splatinsert1230 = insertelement <2 x double> poison, double %251, i32 0
      %broadcast.splat1231 = shufflevector <2 x double> %broadcast.splatinsert1230, <2 x double> poison, <2 x i32> zeroinitializer
      %252 = getelementptr inbounds double, double* %241, i64 %index1226
      %253 = bitcast double* %252 to <2 x double>*
      store <2 x double> %broadcast.splat1231, <2 x double>* %253, align 8
      %254 = getelementptr inbounds double, double* %252, i64 2
      %255 = bitcast double* %254 to <2 x double>*
      store <2 x double> %broadcast.splat1231, <2 x double>* %255, align 8
      %index.next1227 = add nuw i64 %index1226, 4
      %256 = icmp eq i64 %index.next1227, %n.vec1225
      br i1 %256, label %middle.block1219, label %vector.body1221

middle.block1219:                                 ; preds = %vector.body1221
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1229 = icmp eq i64 %214, %n.vec1225
     br i1 %cmp.n1229, label %L639, label %L595.us

L595.us:                                          ; preds = %L595.us, %middle.block1219, %vector.memcheck1210, %L595.us.preheader
     %value_phi169473.us = phi i64 [ %259, %L595.us ], [ %n.vec1225, %middle.block1219 ], [ 0, %L595.us.preheader ], [ 0, %vector.memcheck1210 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %257 = load double, double* %239, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %258 = getelementptr inbounds double, double* %241, i64 %value_phi169473.us
      store double %257, double* %258, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %259 = add nuw nsw i64 %value_phi169473.us, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond846.not = icmp eq i64 %259, %214
; │││└
     br i1 %exitcond846.not, label %L639, label %L595.us

L595:                                             ; preds = %L595, %middle.block1197, %vector.memcheck1188, %L595.preheader
     %value_phi169473 = phi i64 [ %263, %L595 ], [ %n.vec1203, %middle.block1197 ], [ 0, %L595.preheader ], [ 0, %vector.memcheck1188 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %260 = getelementptr inbounds double, double* %239, i64 %value_phi169473
          %261 = load double, double* %260, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %262 = getelementptr inbounds double, double* %241, i64 %value_phi169473
      store double %261, double* %262, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %263 = add nuw nsw i64 %value_phi169473, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond840.not = icmp eq i64 %263, %214
; │││└
     br i1 %exitcond840.not, label %L639, label %L595

L639:                                             ; preds = %L595, %L595.us, %middle.block1219, %middle.block1197, %L574, %L541
; └└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:18 within `advance_check2`
; ┌ @ Base.jl:41 within `dotgetproperty`
; │┌ @ Base.jl:38 within `getproperty`
    %264 = getelementptr inbounds { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }, { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }* %1, i64 0, i32 8
    %265 = load atomic {}*, {}** %264 unordered, align 8
; └└
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ abstractarray.jl:95 within `axes`
; ││┌ @ array.jl:151 within `size`
     %266 = bitcast {}* %265 to { i8*, i64, i16, i16, i32 }*
     %267 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %266, i64 0, i32 1
     %268 = load i64, i64* %267, align 8
; │└└
; │┌ @ broadcast.jl:284 within `instantiate`
; ││┌ @ broadcast.jl:543 within `check_broadcast_axes`
; │││┌ @ abstractarray.jl:95 within `axes`
; ││││┌ @ array.jl:151 within `size`
       %269 = load i64, i64* %104, align 8
; │││└└
; │││┌ @ broadcast.jl:540 within `check_broadcast_shape`
; ││││┌ @ broadcast.jl:518 within `_bcsm`
; │││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
        %270 = icmp ne i64 %268, %269
; │││││└
       %271 = icmp ne i64 %269, 1
; ││││└
      %272 = and i1 %270, %271
      br i1 %272, label %L658, label %L677

L658:                                             ; preds = %L639
      %ptls_field16091647 = getelementptr inbounds {}**, {}*** %9, i64 2
      %273 = bitcast {}*** %ptls_field16091647 to i8**
      %ptls_load161016481649 = load i8*, i8** %273, align 8
      %274 = call noalias nonnull {}* @ijl_gc_pool_alloc(i8* %ptls_load161016481649, i32 1392, i32 16) #7
      %275 = bitcast {}* %274 to i64*
      %276 = getelementptr inbounds i64, i64* %275, i64 -1
      store atomic i64 4712015424, i64* %276 unordered, align 8
      %277 = bitcast {}* %274 to {}**
      store {}* inttoptr (i64 4725386192 to {}*), {}** %277, align 8
      call void @ijl_throw({}* %274)
      unreachable

L677:                                             ; preds = %L639
; │└└└
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:953
; ││┌ @ tuple.jl:398 within `==`
; │││┌ @ tuple.jl:402 within `_eq`
; ││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
       %.not280.not = icmp eq i64 %268, %269
; ││└└└
    br i1 %.not280.not, label %L693, label %L696

L693:                                             ; preds = %L677
    %278 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 9
    store {}* %41, {}** %278, align 8
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954 @ array.jl:346 @ array.jl:322
    %279 = call nonnull {}* @"j__copyto_impl!_2819"({}* nonnull %265, i64 signext 1, {}* nonnull %36, i64 signext 1, i64 signext %268) #0
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954
    br label %L791

L696:                                             ; preds = %L677
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:934 within `broadcast_unalias`
        %.not351 = icmp eq {}* %265, %36
        br i1 %.not351, label %L726, label %L699

L699:                                             ; preds = %L696
; ││││││┌ @ abstractarray.jl:1427 within `unalias`
; │││││││┌ @ abstractarray.jl:1462 within `mightalias`
; ││││││││┌ @ reflection.jl:593 within `isbits`
; │││││││││┌ @ Base.jl:33 within `getproperty`
            %280 = load i8, i8* inttoptr (i64 4708472424 to i8*), align 8
; ││││││││└└
          %281 = and i8 %280, 8
          %.not355.not = icmp eq i8 %281, 0
          br i1 %.not355.not, label %L709, label %L726

L709:                                             ; preds = %L699
; ││││││││┌ @ abstractarray.jl:1486 within `dataids`
; │││││││││┌ @ abstractarray.jl:1187 within `pointer`
; ││││││││││┌ @ pointer.jl:65 within `unsafe_convert`
             %282 = bitcast {}* %265 to i8**
             %283 = load i8*, i8** %282, align 8
             %284 = bitcast {}* %36 to i8**
             %285 = load i8*, i8** %284, align 8
; ││││││││└└└
; ││││││││┌ @ abstractarray.jl:1469 within `_isdisjoint`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %286 = icmp eq i8* %283, %285
; │││││││└└└└
         br i1 %286, label %L721, label %L726

L721:                                             ; preds = %L709
         %287 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 9
         store {}* %41, {}** %287, align 8
; │││││││┌ @ abstractarray.jl:1443 within `unaliascopy`
; ││││││││┌ @ array.jl:369 within `copy`
           %288 = call nonnull {}* inttoptr (i64 4320602180 to {}* ({}*)*)({}* nonnull %36)
; └└└└└└└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:39 within `advance_check2`
; ┌ @ reducedim.jl:1007 within `all`
; │┌ @ reducedim.jl:1007 within `#all#795`
; ││┌ @ reduce.jl:1260 within `_all`
; │││┌ @ array.jl:898 within `iterate`
      br label %L726

L726:                                             ; preds = %L721, %L709, %L699, %L696
      %value_phi158 = phi {}* [ %36, %L696 ], [ %288, %L721 ], [ %36, %L709 ], [ %36, %L699 ]
; └└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:18 within `advance_check2`
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:72 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %.not352 = icmp eq i64 %268, 0
; │││└
     br i1 %.not352, label %L791, label %L747.lr.ph

L747.lr.ph:                                       ; preds = %L726
; ││└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:637 within `extrude`
; ││││││┌ @ broadcast.jl:586 within `newindexer`
; │││││││┌ @ abstractarray.jl:95 within `axes`
; ││││││││┌ @ array.jl:151 within `size`
           %289 = bitcast {}* %value_phi158 to { i8*, i64, i16, i16, i32 }*
           %290 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %289, i64 0, i32 1
           %291 = load i64, i64* %290, align 8
; │││││││└└
; │││││││┌ @ broadcast.jl:587 within `shapeindexer`
; ││││││││┌ @ broadcast.jl:592 within `_newindexer`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %.not354 = icmp eq i64 %291, 1
             %292 = bitcast {}* %value_phi158 to double**
             %293 = load double*, double** %292, align 8
             %294 = bitcast {}* %265 to double**
             %295 = load double*, double** %294, align 8
; ││└└└└└└└└└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:75 within `macro expansion`
     %min.iters.check1268 = icmp ult i64 %268, 4
     br i1 %.not354, label %L747.us.preheader, label %L747.preheader

L747.preheader:                                   ; preds = %L747.lr.ph
     br i1 %min.iters.check1268, label %L747, label %vector.memcheck1234

vector.memcheck1234:                              ; preds = %L747.preheader
     %scevgep1235 = getelementptr double, double* %295, i64 %268
     %scevgep1237 = getelementptr double, double* %293, i64 %268
     %bound01239 = icmp ult double* %295, %scevgep1237
     %bound11240 = icmp ult double* %293, %scevgep1235
     %found.conflict1241 = and i1 %bound01239, %bound11240
     br i1 %found.conflict1241, label %L747, label %vector.ph1247

vector.ph1247:                                    ; preds = %vector.memcheck1234
     %n.vec1249 = and i64 %268, 9223372036854775804
     br label %vector.body1245

vector.body1245:                                  ; preds = %vector.body1245, %vector.ph1247
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1250 = phi i64 [ 0, %vector.ph1247 ], [ %index.next1251, %vector.body1245 ]
      %296 = getelementptr inbounds double, double* %293, i64 %index1250
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %297 = bitcast double* %296 to <2 x double>*
          %wide.load1254 = load <2 x double>, <2 x double>* %297, align 8
          %298 = getelementptr inbounds double, double* %296, i64 2
          %299 = bitcast double* %298 to <2 x double>*
          %wide.load1255 = load <2 x double>, <2 x double>* %299, align 8
; │││└└└└└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %300 = getelementptr inbounds double, double* %295, i64 %index1250
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ array.jl:966 within `setindex!`
      %301 = bitcast double* %300 to <2 x double>*
      store <2 x double> %wide.load1254, <2 x double>* %301, align 8
      %302 = getelementptr inbounds double, double* %300, i64 2
      %303 = bitcast double* %302 to <2 x double>*
      store <2 x double> %wide.load1255, <2 x double>* %303, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index.next1251 = add nuw i64 %index1250, 4
      %304 = icmp eq i64 %index.next1251, %n.vec1249
      br i1 %304, label %middle.block1243, label %vector.body1245

middle.block1243:                                 ; preds = %vector.body1245
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1253 = icmp eq i64 %268, %n.vec1249
     br i1 %cmp.n1253, label %L791, label %L747

L747.us.preheader:                                ; preds = %L747.lr.ph
     br i1 %min.iters.check1268, label %L747.us, label %vector.memcheck1256

vector.memcheck1256:                              ; preds = %L747.us.preheader
     %scevgep1257 = getelementptr double, double* %295, i64 %268
     %scevgep1259 = getelementptr double, double* %293, i64 1
     %bound01261 = icmp ult double* %295, %scevgep1259
     %bound11262 = icmp ult double* %293, %scevgep1257
     %found.conflict1263 = and i1 %bound01261, %bound11262
     br i1 %found.conflict1263, label %L747.us, label %vector.ph1269

vector.ph1269:                                    ; preds = %vector.memcheck1256
     %n.vec1271 = and i64 %268, 9223372036854775804
     br label %vector.body1267

vector.body1267:                                  ; preds = %vector.body1267, %vector.ph1269
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1272 = phi i64 [ 0, %vector.ph1269 ], [ %index.next1273, %vector.body1267 ]
      %305 = load double, double* %293, align 8
      %broadcast.splatinsert1276 = insertelement <2 x double> poison, double %305, i32 0
      %broadcast.splat1277 = shufflevector <2 x double> %broadcast.splatinsert1276, <2 x double> poison, <2 x i32> zeroinitializer
      %306 = getelementptr inbounds double, double* %295, i64 %index1272
      %307 = bitcast double* %306 to <2 x double>*
      store <2 x double> %broadcast.splat1277, <2 x double>* %307, align 8
      %308 = getelementptr inbounds double, double* %306, i64 2
      %309 = bitcast double* %308 to <2 x double>*
      store <2 x double> %broadcast.splat1277, <2 x double>* %309, align 8
      %index.next1273 = add nuw i64 %index1272, 4
      %310 = icmp eq i64 %index.next1273, %n.vec1271
      br i1 %310, label %middle.block1265, label %vector.body1267

middle.block1265:                                 ; preds = %vector.body1267
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1275 = icmp eq i64 %268, %n.vec1271
     br i1 %cmp.n1275, label %L791, label %L747.us

L747.us:                                          ; preds = %L747.us, %middle.block1265, %vector.memcheck1256, %L747.us.preheader
     %value_phi159471.us = phi i64 [ %313, %L747.us ], [ %n.vec1271, %middle.block1265 ], [ 0, %L747.us.preheader ], [ 0, %vector.memcheck1256 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %311 = load double, double* %293, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %312 = getelementptr inbounds double, double* %295, i64 %value_phi159471.us
      store double %311, double* %312, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %313 = add nuw nsw i64 %value_phi159471.us, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond845.not = icmp eq i64 %313, %268
; │││└
     br i1 %exitcond845.not, label %L791, label %L747.us

L747:                                             ; preds = %L747, %middle.block1243, %vector.memcheck1234, %L747.preheader
     %value_phi159471 = phi i64 [ %317, %L747 ], [ %n.vec1249, %middle.block1243 ], [ 0, %L747.preheader ], [ 0, %vector.memcheck1234 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %314 = getelementptr inbounds double, double* %293, i64 %value_phi159471
          %315 = load double, double* %314, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %316 = getelementptr inbounds double, double* %295, i64 %value_phi159471
      store double %315, double* %316, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %317 = add nuw nsw i64 %value_phi159471, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond839.not = icmp eq i64 %317, %268
; │││└
     br i1 %exitcond839.not, label %L791, label %L747

L791:                                             ; preds = %L747, %L747.us, %middle.block1265, %middle.block1243, %L726, %L693
; └└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:19 within `advance_check2`
; ┌ @ Base.jl:41 within `dotgetproperty`
; │┌ @ Base.jl:38 within `getproperty`
    %318 = getelementptr inbounds { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }, { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }* %1, i64 0, i32 9
    %319 = load atomic {}*, {}** %318 unordered, align 8
; └└
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ abstractarray.jl:95 within `axes`
; ││┌ @ array.jl:151 within `size`
     %320 = bitcast {}* %319 to { i8*, i64, i16, i16, i32 }*
     %321 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %320, i64 0, i32 1
     %322 = load i64, i64* %321, align 8
; │└└
; │┌ @ broadcast.jl:284 within `instantiate`
; ││┌ @ broadcast.jl:543 within `check_broadcast_axes`
; │││┌ @ abstractarray.jl:95 within `axes`
; ││││┌ @ array.jl:151 within `size`
       %323 = load i64, i64* %160, align 8
; │││└└
; │││┌ @ broadcast.jl:540 within `check_broadcast_shape`
; ││││┌ @ broadcast.jl:518 within `_bcsm`
; │││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
        %324 = icmp ne i64 %322, %323
; │││││└
       %325 = icmp ne i64 %323, 1
; ││││└
      %326 = and i1 %324, %325
      br i1 %326, label %L810, label %L829

L810:                                             ; preds = %L791
      %ptls_field16111644 = getelementptr inbounds {}**, {}*** %9, i64 2
      %327 = bitcast {}*** %ptls_field16111644 to i8**
      %ptls_load161216451646 = load i8*, i8** %327, align 8
      %328 = call noalias nonnull {}* @ijl_gc_pool_alloc(i8* %ptls_load161216451646, i32 1392, i32 16) #7
      %329 = bitcast {}* %328 to i64*
      %330 = getelementptr inbounds i64, i64* %329, i64 -1
      store atomic i64 4712015424, i64* %330 unordered, align 8
      %331 = bitcast {}* %328 to {}**
      store {}* inttoptr (i64 4725386192 to {}*), {}** %331, align 8
      call void @ijl_throw({}* %328)
      unreachable

L829:                                             ; preds = %L791
; │└└└
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:953
; ││┌ @ tuple.jl:398 within `==`
; │││┌ @ tuple.jl:402 within `_eq`
; ││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
       %.not282.not = icmp eq i64 %322, %323
; ││└└└
    br i1 %.not282.not, label %L845, label %L848

L845:                                             ; preds = %L829
    %332 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 9
    store {}* %41, {}** %332, align 8
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954 @ array.jl:346 @ array.jl:322
    %333 = call nonnull {}* @"j__copyto_impl!_2820"({}* nonnull %319, i64 signext 1, {}* nonnull %41, i64 signext 1, i64 signext %322) #0
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954
    br label %L943.L944_crit_edge

L848:                                             ; preds = %L829
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:934 within `broadcast_unalias`
        %.not344 = icmp eq {}* %319, %41
        br i1 %.not344, label %L878, label %L851

L851:                                             ; preds = %L848
; ││││││┌ @ abstractarray.jl:1427 within `unalias`
; │││││││┌ @ abstractarray.jl:1462 within `mightalias`
; ││││││││┌ @ reflection.jl:593 within `isbits`
; │││││││││┌ @ Base.jl:33 within `getproperty`
            %334 = load i8, i8* inttoptr (i64 4708472424 to i8*), align 8
; ││││││││└└
          %335 = and i8 %334, 8
          %.not348.not = icmp eq i8 %335, 0
          br i1 %.not348.not, label %L861, label %L878

L861:                                             ; preds = %L851
; ││││││││┌ @ abstractarray.jl:1486 within `dataids`
; │││││││││┌ @ abstractarray.jl:1187 within `pointer`
; ││││││││││┌ @ pointer.jl:65 within `unsafe_convert`
             %336 = bitcast {}* %319 to i8**
             %337 = load i8*, i8** %336, align 8
             %338 = bitcast {}* %41 to i8**
             %339 = load i8*, i8** %338, align 8
; ││││││││└└└
; ││││││││┌ @ abstractarray.jl:1469 within `_isdisjoint`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %340 = icmp eq i8* %337, %339
; │││││││└└└└
         br i1 %340, label %L873, label %L878

L873:                                             ; preds = %L861
         %341 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 9
         store {}* %41, {}** %341, align 8
; │││││││┌ @ abstractarray.jl:1443 within `unaliascopy`
; ││││││││┌ @ array.jl:369 within `copy`
           %342 = call nonnull {}* inttoptr (i64 4320602180 to {}* ({}*)*)({}* nonnull %41)
; └└└└└└└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:39 within `advance_check2`
; ┌ @ reducedim.jl:1007 within `all`
; │┌ @ reducedim.jl:1007 within `#all#795`
; ││┌ @ reduce.jl:1260 within `_all`
; │││┌ @ array.jl:898 within `iterate`
      br label %L878

L878:                                             ; preds = %L873, %L861, %L851, %L848
      %value_phi148 = phi {}* [ %41, %L848 ], [ %342, %L873 ], [ %41, %L861 ], [ %41, %L851 ]
; └└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:19 within `advance_check2`
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:72 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %.not345 = icmp eq i64 %322, 0
; │││└
     br i1 %.not345, label %L943.L944_crit_edge, label %L899.lr.ph

L899.lr.ph:                                       ; preds = %L878
; ││└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:637 within `extrude`
; ││││││┌ @ broadcast.jl:586 within `newindexer`
; │││││││┌ @ abstractarray.jl:95 within `axes`
; ││││││││┌ @ array.jl:151 within `size`
           %343 = bitcast {}* %value_phi148 to { i8*, i64, i16, i16, i32 }*
           %344 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %343, i64 0, i32 1
           %345 = load i64, i64* %344, align 8
; │││││││└└
; │││││││┌ @ broadcast.jl:587 within `shapeindexer`
; ││││││││┌ @ broadcast.jl:592 within `_newindexer`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %.not347 = icmp eq i64 %345, 1
             %346 = bitcast {}* %value_phi148 to double**
             %347 = load double*, double** %346, align 8
             %348 = bitcast {}* %319 to double**
             %349 = load double*, double** %348, align 8
; ││└└└└└└└└└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:75 within `macro expansion`
     %min.iters.check1314 = icmp ult i64 %322, 4
     br i1 %.not347, label %L899.us.preheader, label %L899.preheader

L899.preheader:                                   ; preds = %L899.lr.ph
     br i1 %min.iters.check1314, label %L899, label %vector.memcheck1280

vector.memcheck1280:                              ; preds = %L899.preheader
     %scevgep1281 = getelementptr double, double* %349, i64 %322
     %scevgep1283 = getelementptr double, double* %347, i64 %322
     %bound01285 = icmp ult double* %349, %scevgep1283
     %bound11286 = icmp ult double* %347, %scevgep1281
     %found.conflict1287 = and i1 %bound01285, %bound11286
     br i1 %found.conflict1287, label %L899, label %vector.ph1293

vector.ph1293:                                    ; preds = %vector.memcheck1280
     %n.vec1295 = and i64 %322, 9223372036854775804
     br label %vector.body1291

vector.body1291:                                  ; preds = %vector.body1291, %vector.ph1293
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1296 = phi i64 [ 0, %vector.ph1293 ], [ %index.next1297, %vector.body1291 ]
      %350 = getelementptr inbounds double, double* %347, i64 %index1296
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %351 = bitcast double* %350 to <2 x double>*
          %wide.load1300 = load <2 x double>, <2 x double>* %351, align 8
          %352 = getelementptr inbounds double, double* %350, i64 2
          %353 = bitcast double* %352 to <2 x double>*
          %wide.load1301 = load <2 x double>, <2 x double>* %353, align 8
; │││└└└└└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %354 = getelementptr inbounds double, double* %349, i64 %index1296
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ array.jl:966 within `setindex!`
      %355 = bitcast double* %354 to <2 x double>*
      store <2 x double> %wide.load1300, <2 x double>* %355, align 8
      %356 = getelementptr inbounds double, double* %354, i64 2
      %357 = bitcast double* %356 to <2 x double>*
      store <2 x double> %wide.load1301, <2 x double>* %357, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index.next1297 = add nuw i64 %index1296, 4
      %358 = icmp eq i64 %index.next1297, %n.vec1295
      br i1 %358, label %middle.block1289, label %vector.body1291

middle.block1289:                                 ; preds = %vector.body1291
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1299 = icmp eq i64 %322, %n.vec1295
     br i1 %cmp.n1299, label %L943.L944_crit_edge, label %L899

L899.us.preheader:                                ; preds = %L899.lr.ph
     br i1 %min.iters.check1314, label %L899.us, label %vector.memcheck1302

vector.memcheck1302:                              ; preds = %L899.us.preheader
     %scevgep1303 = getelementptr double, double* %349, i64 %322
     %scevgep1305 = getelementptr double, double* %347, i64 1
     %bound01307 = icmp ult double* %349, %scevgep1305
     %bound11308 = icmp ult double* %347, %scevgep1303
     %found.conflict1309 = and i1 %bound01307, %bound11308
     br i1 %found.conflict1309, label %L899.us, label %vector.ph1315

vector.ph1315:                                    ; preds = %vector.memcheck1302
     %n.vec1317 = and i64 %322, 9223372036854775804
     br label %vector.body1313

vector.body1313:                                  ; preds = %vector.body1313, %vector.ph1315
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1318 = phi i64 [ 0, %vector.ph1315 ], [ %index.next1319, %vector.body1313 ]
      %359 = load double, double* %347, align 8
      %broadcast.splatinsert1322 = insertelement <2 x double> poison, double %359, i32 0
      %broadcast.splat1323 = shufflevector <2 x double> %broadcast.splatinsert1322, <2 x double> poison, <2 x i32> zeroinitializer
      %360 = getelementptr inbounds double, double* %349, i64 %index1318
      %361 = bitcast double* %360 to <2 x double>*
      store <2 x double> %broadcast.splat1323, <2 x double>* %361, align 8
      %362 = getelementptr inbounds double, double* %360, i64 2
      %363 = bitcast double* %362 to <2 x double>*
      store <2 x double> %broadcast.splat1323, <2 x double>* %363, align 8
      %index.next1319 = add nuw i64 %index1318, 4
      %364 = icmp eq i64 %index.next1319, %n.vec1317
      br i1 %364, label %middle.block1311, label %vector.body1313

middle.block1311:                                 ; preds = %vector.body1313
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1321 = icmp eq i64 %322, %n.vec1317
     br i1 %cmp.n1321, label %L943.L944_crit_edge, label %L899.us

L899.us:                                          ; preds = %L899.us, %middle.block1311, %vector.memcheck1302, %L899.us.preheader
     %value_phi149469.us = phi i64 [ %367, %L899.us ], [ %n.vec1317, %middle.block1311 ], [ 0, %L899.us.preheader ], [ 0, %vector.memcheck1302 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %365 = load double, double* %347, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %366 = getelementptr inbounds double, double* %349, i64 %value_phi149469.us
      store double %365, double* %366, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %367 = add nuw nsw i64 %value_phi149469.us, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond844.not = icmp eq i64 %367, %322
; │││└
     br i1 %exitcond844.not, label %L943.L944_crit_edge, label %L899.us

L899:                                             ; preds = %L899, %middle.block1289, %vector.memcheck1280, %L899.preheader
     %value_phi149469 = phi i64 [ %371, %L899 ], [ %n.vec1295, %middle.block1289 ], [ 0, %L899.preheader ], [ 0, %vector.memcheck1280 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %368 = getelementptr inbounds double, double* %347, i64 %value_phi149469
          %369 = load double, double* %368, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %370 = getelementptr inbounds double, double* %349, i64 %value_phi149469
      store double %369, double* %370, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %371 = add nuw nsw i64 %value_phi149469, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond838.not = icmp eq i64 %371, %322
; │││└
     br i1 %exitcond838.not, label %L943.L944_crit_edge, label %L899

L943.L944_crit_edge:                              ; preds = %L899, %L899.us, %middle.block1311, %middle.block1289, %L878, %L845
     %.pre-phi.pre-phi.pre-phi = phi {}* [ %319, %middle.block1289 ], [ %319, %L899 ], [ %319, %middle.block1311 ], [ %319, %L899.us ], [ %319, %L878 ], [ %319, %L845 ]
     %372 = bitcast {}* %25 to { i8*, i64, i16, i16, i32 }*
     %373 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %372, i64 0, i32 1
     %374 = bitcast {}* %25 to double**
     %375 = getelementptr inbounds { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }, { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }* %1, i64 0, i32 33
     %376 = load atomic {}*, {}** %375 unordered, align 8
     %377 = load double, double* %15, align 8
     %378 = bitcast {}* %211 to double**
     %379 = getelementptr inbounds { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }, { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }* %1, i64 0, i32 35
     %380 = load atomic {}*, {}** %379 unordered, align 8
     %381 = bitcast {}* %265 to double**
     %382 = getelementptr inbounds { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }, { i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }* %1, i64 0, i32 38
     %383 = load atomic {}*, {}** %382 unordered, align 8
     %384 = bitcast {}* %319 to double**
     %385 = bitcast {}* %16 to { i8*, i64, i16, i16, i32 }*
     %386 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %385, i64 0, i32 1
     %387 = bitcast {}* %43 to double**
     %388 = bitcast {}* %99 to double**
     %389 = bitcast {}* %155 to double**
; └└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:6 within `advance_check2`
; ┌ @ array.jl:126 within `vect`
; │┌ @ range.jl:883 within `iterate`
    br label %L944

L944:                                             ; preds = %L1895, %L943.L944_crit_edge
    %value_phi24 = phi i64 [ 1, %L943.L944_crit_edge ], [ %659, %L1895 ]
    %390 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 9
    store {}* %41, {}** %390, align 8
    %391 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 5
    store {}* %383, {}** %391, align 8
    %392 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 4
    store {}* %380, {}** %392, align 16
    %393 = getelementptr inbounds [11 x {}*], [11 x {}*]* %gcframe1625, i64 0, i64 3
    store {}* %376, {}** %393, align 8
; └└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:23 within `advance_check2`
  call void @j_comp_u_v_eta_t_check2_2821({ i64, i64, i64, i64, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}*, {}* }* nocapture nonnull readonly %1, { double, double, double, double, double, double, double, double, {}*, {}* }* nocapture nonnull readonly %2, { i64, i64, i64, i64, i64, i64, i64, i64, double, double }* nocapture nonnull readonly %3, [12 x { i64, i64, {}*, {}*, {}* }]* nocapture nonnull readonly %4, [14 x { i64, i64, {}*, {}*, {}* }]* nocapture nonnull readonly %5, { { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, {}*, {}*, {}*, {}*, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* }, { i64, i64, {}*, {}*, {}* } }* nocapture nonnull readonly %6) #0
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:25 within `advance_check2`
; ┌ @ int.jl:83 within `<`
   %394 = icmp ugt i64 %value_phi24, 3
; └
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl within `advance_check2`
  %.pre850 = add nsw i64 %value_phi24, -1
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:25 within `advance_check2`
  br i1 %394, label %L1421, label %L949

L949:                                             ; preds = %L944
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:26 within `advance_check2`
; ┌ @ array.jl:924 within `getindex`
   %395 = load i64, i64* %373, align 8
   %396 = icmp ult i64 %.pre850, %395
   br i1 %396, label %idxend, label %oob

L973:                                             ; preds = %idxend
; └
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ broadcast.jl:284 within `instantiate`
; ││┌ @ broadcast.jl:543 within `check_broadcast_axes`
; │││┌ @ broadcast.jl:540 within `check_broadcast_shape`
      %ptls_field16131641 = getelementptr inbounds {}**, {}*** %9, i64 2
      %397 = bitcast {}*** %ptls_field16131641 to i8**
      %ptls_load161416421643 = load i8*, i8** %397, align 8
      %398 = call noalias nonnull {}* @ijl_gc_pool_alloc(i8* %ptls_load161416421643, i32 1392, i32 16) #7
      %399 = bitcast {}* %398 to i64*
      %400 = getelementptr inbounds i64, i64* %399, i64 -1
      store atomic i64 4712015424, i64* %400 unordered, align 8
      %401 = bitcast {}* %398 to {}**
      store {}* inttoptr (i64 4725386192 to {}*), {}** %401, align 8
      call void @ijl_throw({}* %398)
      unreachable

L992:                                             ; preds = %idxend
; │└└└
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:953
; ││┌ @ tuple.jl:398 within `==`
; │││┌ @ tuple.jl:402 within `_eq`
; ││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
       %.not284.not = icmp eq i64 %710, %713
; ││└└└
    br i1 %.not284.not, label %L1008, label %L1011

L1008:                                            ; preds = %L992
    store {}* %709, {}** %28, align 16
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954 @ array.jl:346 @ array.jl:322
    %402 = call nonnull {}* @"j__copyto_impl!_2824"({}* nonnull %211, i64 signext 1, {}* nonnull %709, i64 signext 1, i64 signext %710) #0
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954
    br label %L1106

L1011:                                            ; preds = %L992
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:934 within `broadcast_unalias`
        %.not337 = icmp eq {}* %211, %709
        br i1 %.not337, label %L1041, label %L1014

L1014:                                            ; preds = %L1011
; ││││││┌ @ abstractarray.jl:1427 within `unalias`
; │││││││┌ @ abstractarray.jl:1462 within `mightalias`
; ││││││││┌ @ reflection.jl:593 within `isbits`
; │││││││││┌ @ Base.jl:33 within `getproperty`
            %403 = load i8, i8* inttoptr (i64 4708472424 to i8*), align 8
; ││││││││└└
          %404 = and i8 %403, 8
          %.not341.not = icmp eq i8 %404, 0
          br i1 %.not341.not, label %L1024, label %L1041

L1024:                                            ; preds = %L1014
; ││││││││┌ @ abstractarray.jl:1486 within `dataids`
; │││││││││┌ @ abstractarray.jl:1187 within `pointer`
; ││││││││││┌ @ pointer.jl:65 within `unsafe_convert`
             %405 = bitcast {}* %211 to i8**
             %406 = load i8*, i8** %405, align 8
             %407 = bitcast {}* %709 to i8**
             %408 = load i8*, i8** %407, align 8
; ││││││││└└└
; ││││││││┌ @ abstractarray.jl:1469 within `_isdisjoint`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %409 = icmp eq i8* %406, %408
; │││││││└└└└
         br i1 %409, label %L1036, label %L1041

L1036:                                            ; preds = %L1024
         store {}* %709, {}** %28, align 16
; │││││││┌ @ abstractarray.jl:1443 within `unaliascopy`
; ││││││││┌ @ array.jl:369 within `copy`
           %410 = call nonnull {}* inttoptr (i64 4320602180 to {}* ({}*)*)({}* nonnull %709)
; └└└└└└└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:39 within `advance_check2`
; ┌ @ reducedim.jl:1007 within `all`
; │┌ @ reducedim.jl:1007 within `#all#795`
; ││┌ @ reduce.jl:1260 within `_all`
; │││┌ @ array.jl:898 within `iterate`
      br label %L1041

L1041:                                            ; preds = %L1036, %L1024, %L1014, %L1011
      %value_phi138 = phi {}* [ %709, %L1011 ], [ %410, %L1036 ], [ %709, %L1024 ], [ %709, %L1014 ]
; └└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:26 within `advance_check2`
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:72 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %.not338 = icmp eq i64 %710, 0
; │││└
     br i1 %.not338, label %L1106, label %L1062.lr.ph

L1062.lr.ph:                                      ; preds = %L1041
; ││└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:637 within `extrude`
; ││││││┌ @ broadcast.jl:586 within `newindexer`
; │││││││┌ @ abstractarray.jl:95 within `axes`
; ││││││││┌ @ array.jl:151 within `size`
           %411 = bitcast {}* %value_phi138 to { i8*, i64, i16, i16, i32 }*
           %412 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %411, i64 0, i32 1
           %413 = load i64, i64* %412, align 8
; │││││││└└
; │││││││┌ @ broadcast.jl:587 within `shapeindexer`
; ││││││││┌ @ broadcast.jl:592 within `_newindexer`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %.not340 = icmp eq i64 %413, 1
             %414 = bitcast {}* %value_phi138 to double**
             %415 = load double*, double** %414, align 8
             %416 = load double*, double** %378, align 8
; ││└└└└└└└└└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:75 within `macro expansion`
     %min.iters.check1568 = icmp ult i64 %710, 4
     br i1 %.not340, label %L1062.us.preheader, label %L1062.preheader

L1062.preheader:                                  ; preds = %L1062.lr.ph
     br i1 %min.iters.check1568, label %L1062, label %vector.memcheck1580

vector.memcheck1580:                              ; preds = %L1062.preheader
     %scevgep1581 = getelementptr double, double* %416, i64 %710
     %scevgep1583 = getelementptr double, double* %415, i64 %710
     %bound01585 = icmp ult double* %416, %scevgep1583
     %bound11586 = icmp ult double* %415, %scevgep1581
     %found.conflict1587 = and i1 %bound01585, %bound11586
     br i1 %found.conflict1587, label %L1062, label %vector.ph1593

vector.ph1593:                                    ; preds = %vector.memcheck1580
     %n.vec1595 = and i64 %710, 9223372036854775804
     br label %vector.body1591

vector.body1591:                                  ; preds = %vector.body1591, %vector.ph1593
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1596 = phi i64 [ 0, %vector.ph1593 ], [ %index.next1597, %vector.body1591 ]
      %417 = getelementptr inbounds double, double* %415, i64 %index1596
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %418 = bitcast double* %417 to <2 x double>*
          %wide.load1600 = load <2 x double>, <2 x double>* %418, align 8
          %419 = getelementptr inbounds double, double* %417, i64 2
          %420 = bitcast double* %419 to <2 x double>*
          %wide.load1601 = load <2 x double>, <2 x double>* %420, align 8
; │││└└└└└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %421 = getelementptr inbounds double, double* %416, i64 %index1596
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ array.jl:966 within `setindex!`
      %422 = bitcast double* %421 to <2 x double>*
      store <2 x double> %wide.load1600, <2 x double>* %422, align 8
      %423 = getelementptr inbounds double, double* %421, i64 2
      %424 = bitcast double* %423 to <2 x double>*
      store <2 x double> %wide.load1601, <2 x double>* %424, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index.next1597 = add nuw i64 %index1596, 4
      %425 = icmp eq i64 %index.next1597, %n.vec1595
      br i1 %425, label %middle.block1589, label %vector.body1591

middle.block1589:                                 ; preds = %vector.body1591
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1599 = icmp eq i64 %710, %n.vec1595
     br i1 %cmp.n1599, label %L1106, label %L1062

L1062.us.preheader:                               ; preds = %L1062.lr.ph
     br i1 %min.iters.check1568, label %L1062.us, label %vector.memcheck1556

vector.memcheck1556:                              ; preds = %L1062.us.preheader
     %scevgep1557 = getelementptr double, double* %416, i64 %710
     %scevgep1559 = getelementptr double, double* %415, i64 1
     %bound01561 = icmp ult double* %416, %scevgep1559
     %bound11562 = icmp ult double* %415, %scevgep1557
     %found.conflict1563 = and i1 %bound01561, %bound11562
     br i1 %found.conflict1563, label %L1062.us, label %vector.ph1569

vector.ph1569:                                    ; preds = %vector.memcheck1556
     %n.vec1571 = and i64 %710, 9223372036854775804
     br label %vector.body1567

vector.body1567:                                  ; preds = %vector.body1567, %vector.ph1569
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1572 = phi i64 [ 0, %vector.ph1569 ], [ %index.next1573, %vector.body1567 ]
      %426 = load double, double* %415, align 8
      %broadcast.splatinsert1576 = insertelement <2 x double> poison, double %426, i32 0
      %broadcast.splat1577 = shufflevector <2 x double> %broadcast.splatinsert1576, <2 x double> poison, <2 x i32> zeroinitializer
      %427 = getelementptr inbounds double, double* %416, i64 %index1572
      %428 = bitcast double* %427 to <2 x double>*
      store <2 x double> %broadcast.splat1577, <2 x double>* %428, align 8
      %429 = getelementptr inbounds double, double* %427, i64 2
      %430 = bitcast double* %429 to <2 x double>*
      store <2 x double> %broadcast.splat1577, <2 x double>* %430, align 8
      %index.next1573 = add nuw i64 %index1572, 4
      %431 = icmp eq i64 %index.next1573, %n.vec1571
      br i1 %431, label %middle.block1565, label %vector.body1567

middle.block1565:                                 ; preds = %vector.body1567
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1575 = icmp eq i64 %710, %n.vec1571
     br i1 %cmp.n1575, label %L1106, label %L1062.us

L1062.us:                                         ; preds = %L1062.us, %middle.block1565, %vector.memcheck1556, %L1062.us.preheader
     %value_phi139457.us = phi i64 [ %434, %L1062.us ], [ %n.vec1571, %middle.block1565 ], [ 0, %L1062.us.preheader ], [ 0, %vector.memcheck1556 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %432 = load double, double* %415, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %433 = getelementptr inbounds double, double* %416, i64 %value_phi139457.us
      store double %432, double* %433, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %434 = add nuw nsw i64 %value_phi139457.us, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond832.not = icmp eq i64 %434, %710
; │││└
     br i1 %exitcond832.not, label %L1106, label %L1062.us

L1062:                                            ; preds = %L1062, %middle.block1589, %vector.memcheck1580, %L1062.preheader
     %value_phi139457 = phi i64 [ %438, %L1062 ], [ %n.vec1595, %middle.block1589 ], [ 0, %L1062.preheader ], [ 0, %vector.memcheck1580 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %435 = getelementptr inbounds double, double* %415, i64 %value_phi139457
          %436 = load double, double* %435, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %437 = getelementptr inbounds double, double* %416, i64 %value_phi139457
      store double %436, double* %437, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %438 = add nuw nsw i64 %value_phi139457, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond826.not = icmp eq i64 %438, %710
; │││└
     br i1 %exitcond826.not, label %L1106, label %L1062

L1106:                                            ; preds = %L1062, %L1062.us, %middle.block1565, %middle.block1589, %L1041, %L1008
; └└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:27 within `advance_check2`
; ┌ @ array.jl:924 within `getindex`
   %439 = load i64, i64* %373, align 8
   %440 = icmp ult i64 %.pre850, %439
   br i1 %440, label %idxend30, label %oob29

L1130:                                            ; preds = %idxend30
; └
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ broadcast.jl:284 within `instantiate`
; ││┌ @ broadcast.jl:543 within `check_broadcast_axes`
; │││┌ @ broadcast.jl:540 within `check_broadcast_shape`
      %ptls_field16151638 = getelementptr inbounds {}**, {}*** %9, i64 2
      %441 = bitcast {}*** %ptls_field16151638 to i8**
      %ptls_load161616391640 = load i8*, i8** %441, align 8
      %442 = call noalias nonnull {}* @ijl_gc_pool_alloc(i8* %ptls_load161616391640, i32 1392, i32 16) #7
      %443 = bitcast {}* %442 to i64*
      %444 = getelementptr inbounds i64, i64* %443, i64 -1
      store atomic i64 4712015424, i64* %444 unordered, align 8
      %445 = bitcast {}* %442 to {}**
      store {}* inttoptr (i64 4725386192 to {}*), {}** %445, align 8
      call void @ijl_throw({}* %442)
      unreachable

L1149:                                            ; preds = %idxend30
; │└└└
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:953
; ││┌ @ tuple.jl:398 within `==`
; │││┌ @ tuple.jl:402 within `_eq`
; ││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
       %.not286.not = icmp eq i64 %724, %727
; ││└└└
    br i1 %.not286.not, label %L1165, label %L1168

L1165:                                            ; preds = %L1149
    store {}* %723, {}** %28, align 16
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954 @ array.jl:346 @ array.jl:322
    %446 = call nonnull {}* @"j__copyto_impl!_2827"({}* nonnull %265, i64 signext 1, {}* nonnull %723, i64 signext 1, i64 signext %724) #0
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954
    br label %L1263

L1168:                                            ; preds = %L1149
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:934 within `broadcast_unalias`
        %.not330 = icmp eq {}* %265, %723
        br i1 %.not330, label %L1198, label %L1171

L1171:                                            ; preds = %L1168
; ││││││┌ @ abstractarray.jl:1427 within `unalias`
; │││││││┌ @ abstractarray.jl:1462 within `mightalias`
; ││││││││┌ @ reflection.jl:593 within `isbits`
; │││││││││┌ @ Base.jl:33 within `getproperty`
            %447 = load i8, i8* inttoptr (i64 4708472424 to i8*), align 8
; ││││││││└└
          %448 = and i8 %447, 8
          %.not334.not = icmp eq i8 %448, 0
          br i1 %.not334.not, label %L1181, label %L1198

L1181:                                            ; preds = %L1171
; ││││││││┌ @ abstractarray.jl:1486 within `dataids`
; │││││││││┌ @ abstractarray.jl:1187 within `pointer`
; ││││││││││┌ @ pointer.jl:65 within `unsafe_convert`
             %449 = bitcast {}* %265 to i8**
             %450 = load i8*, i8** %449, align 8
             %451 = bitcast {}* %723 to i8**
             %452 = load i8*, i8** %451, align 8
; ││││││││└└└
; ││││││││┌ @ abstractarray.jl:1469 within `_isdisjoint`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %453 = icmp eq i8* %450, %452
; │││││││└└└└
         br i1 %453, label %L1193, label %L1198

L1193:                                            ; preds = %L1181
         store {}* %723, {}** %28, align 16
; │││││││┌ @ abstractarray.jl:1443 within `unaliascopy`
; ││││││││┌ @ array.jl:369 within `copy`
           %454 = call nonnull {}* inttoptr (i64 4320602180 to {}* ({}*)*)({}* nonnull %723)
; └└└└└└└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:39 within `advance_check2`
; ┌ @ reducedim.jl:1007 within `all`
; │┌ @ reducedim.jl:1007 within `#all#795`
; ││┌ @ reduce.jl:1260 within `_all`
; │││┌ @ array.jl:898 within `iterate`
      br label %L1198

L1198:                                            ; preds = %L1193, %L1181, %L1171, %L1168
      %value_phi128 = phi {}* [ %723, %L1168 ], [ %454, %L1193 ], [ %723, %L1181 ], [ %723, %L1171 ]
; └└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:27 within `advance_check2`
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:72 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %.not331 = icmp eq i64 %724, 0
; │││└
     br i1 %.not331, label %L1263, label %L1219.lr.ph

L1219.lr.ph:                                      ; preds = %L1198
; ││└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:637 within `extrude`
; ││││││┌ @ broadcast.jl:586 within `newindexer`
; │││││││┌ @ abstractarray.jl:95 within `axes`
; ││││││││┌ @ array.jl:151 within `size`
           %455 = bitcast {}* %value_phi128 to { i8*, i64, i16, i16, i32 }*
           %456 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %455, i64 0, i32 1
           %457 = load i64, i64* %456, align 8
; │││││││└└
; │││││││┌ @ broadcast.jl:587 within `shapeindexer`
; ││││││││┌ @ broadcast.jl:592 within `_newindexer`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %.not333 = icmp eq i64 %457, 1
             %458 = bitcast {}* %value_phi128 to double**
             %459 = load double*, double** %458, align 8
             %460 = load double*, double** %381, align 8
; ││└└└└└└└└└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:75 within `macro expansion`
     %min.iters.check1522 = icmp ult i64 %724, 4
     br i1 %.not333, label %L1219.us.preheader, label %L1219.preheader

L1219.preheader:                                  ; preds = %L1219.lr.ph
     br i1 %min.iters.check1522, label %L1219, label %vector.memcheck1534

vector.memcheck1534:                              ; preds = %L1219.preheader
     %scevgep1535 = getelementptr double, double* %460, i64 %724
     %scevgep1537 = getelementptr double, double* %459, i64 %724
     %bound01539 = icmp ult double* %460, %scevgep1537
     %bound11540 = icmp ult double* %459, %scevgep1535
     %found.conflict1541 = and i1 %bound01539, %bound11540
     br i1 %found.conflict1541, label %L1219, label %vector.ph1547

vector.ph1547:                                    ; preds = %vector.memcheck1534
     %n.vec1549 = and i64 %724, 9223372036854775804
     br label %vector.body1545

vector.body1545:                                  ; preds = %vector.body1545, %vector.ph1547
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1550 = phi i64 [ 0, %vector.ph1547 ], [ %index.next1551, %vector.body1545 ]
      %461 = getelementptr inbounds double, double* %459, i64 %index1550
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %462 = bitcast double* %461 to <2 x double>*
          %wide.load1554 = load <2 x double>, <2 x double>* %462, align 8
          %463 = getelementptr inbounds double, double* %461, i64 2
          %464 = bitcast double* %463 to <2 x double>*
          %wide.load1555 = load <2 x double>, <2 x double>* %464, align 8
; │││└└└└└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %465 = getelementptr inbounds double, double* %460, i64 %index1550
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ array.jl:966 within `setindex!`
      %466 = bitcast double* %465 to <2 x double>*
      store <2 x double> %wide.load1554, <2 x double>* %466, align 8
      %467 = getelementptr inbounds double, double* %465, i64 2
      %468 = bitcast double* %467 to <2 x double>*
      store <2 x double> %wide.load1555, <2 x double>* %468, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index.next1551 = add nuw i64 %index1550, 4
      %469 = icmp eq i64 %index.next1551, %n.vec1549
      br i1 %469, label %middle.block1543, label %vector.body1545

middle.block1543:                                 ; preds = %vector.body1545
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1553 = icmp eq i64 %724, %n.vec1549
     br i1 %cmp.n1553, label %L1263, label %L1219

L1219.us.preheader:                               ; preds = %L1219.lr.ph
     br i1 %min.iters.check1522, label %L1219.us, label %vector.memcheck1510

vector.memcheck1510:                              ; preds = %L1219.us.preheader
     %scevgep1511 = getelementptr double, double* %460, i64 %724
     %scevgep1513 = getelementptr double, double* %459, i64 1
     %bound01515 = icmp ult double* %460, %scevgep1513
     %bound11516 = icmp ult double* %459, %scevgep1511
     %found.conflict1517 = and i1 %bound01515, %bound11516
     br i1 %found.conflict1517, label %L1219.us, label %vector.ph1523

vector.ph1523:                                    ; preds = %vector.memcheck1510
     %n.vec1525 = and i64 %724, 9223372036854775804
     br label %vector.body1521

vector.body1521:                                  ; preds = %vector.body1521, %vector.ph1523
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1526 = phi i64 [ 0, %vector.ph1523 ], [ %index.next1527, %vector.body1521 ]
      %470 = load double, double* %459, align 8
      %broadcast.splatinsert1530 = insertelement <2 x double> poison, double %470, i32 0
      %broadcast.splat1531 = shufflevector <2 x double> %broadcast.splatinsert1530, <2 x double> poison, <2 x i32> zeroinitializer
      %471 = getelementptr inbounds double, double* %460, i64 %index1526
      %472 = bitcast double* %471 to <2 x double>*
      store <2 x double> %broadcast.splat1531, <2 x double>* %472, align 8
      %473 = getelementptr inbounds double, double* %471, i64 2
      %474 = bitcast double* %473 to <2 x double>*
      store <2 x double> %broadcast.splat1531, <2 x double>* %474, align 8
      %index.next1527 = add nuw i64 %index1526, 4
      %475 = icmp eq i64 %index.next1527, %n.vec1525
      br i1 %475, label %middle.block1519, label %vector.body1521

middle.block1519:                                 ; preds = %vector.body1521
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1529 = icmp eq i64 %724, %n.vec1525
     br i1 %cmp.n1529, label %L1263, label %L1219.us

L1219.us:                                         ; preds = %L1219.us, %middle.block1519, %vector.memcheck1510, %L1219.us.preheader
     %value_phi129459.us = phi i64 [ %478, %L1219.us ], [ %n.vec1525, %middle.block1519 ], [ 0, %L1219.us.preheader ], [ 0, %vector.memcheck1510 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %476 = load double, double* %459, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %477 = getelementptr inbounds double, double* %460, i64 %value_phi129459.us
      store double %476, double* %477, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %478 = add nuw nsw i64 %value_phi129459.us, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond833.not = icmp eq i64 %478, %724
; │││└
     br i1 %exitcond833.not, label %L1263, label %L1219.us

L1219:                                            ; preds = %L1219, %middle.block1543, %vector.memcheck1534, %L1219.preheader
     %value_phi129459 = phi i64 [ %482, %L1219 ], [ %n.vec1549, %middle.block1543 ], [ 0, %L1219.preheader ], [ 0, %vector.memcheck1534 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %479 = getelementptr inbounds double, double* %459, i64 %value_phi129459
          %480 = load double, double* %479, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %481 = getelementptr inbounds double, double* %460, i64 %value_phi129459
      store double %480, double* %481, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %482 = add nuw nsw i64 %value_phi129459, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond827.not = icmp eq i64 %482, %724
; │││└
     br i1 %exitcond827.not, label %L1263, label %L1219

L1263:                                            ; preds = %L1219, %L1219.us, %middle.block1519, %middle.block1543, %L1198, %L1165
; └└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:28 within `advance_check2`
; ┌ @ array.jl:924 within `getindex`
   %483 = load i64, i64* %373, align 8
   %484 = icmp ult i64 %.pre850, %483
   br i1 %484, label %idxend35, label %oob34

L1287:                                            ; preds = %idxend35
; └
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ broadcast.jl:284 within `instantiate`
; ││┌ @ broadcast.jl:543 within `check_broadcast_axes`
; │││┌ @ broadcast.jl:540 within `check_broadcast_shape`
      %ptls_field16171635 = getelementptr inbounds {}**, {}*** %9, i64 2
      %485 = bitcast {}*** %ptls_field16171635 to i8**
      %ptls_load161816361637 = load i8*, i8** %485, align 8
      %486 = call noalias nonnull {}* @ijl_gc_pool_alloc(i8* %ptls_load161816361637, i32 1392, i32 16) #7
      %487 = bitcast {}* %486 to i64*
      %488 = getelementptr inbounds i64, i64* %487, i64 -1
      store atomic i64 4712015424, i64* %488 unordered, align 8
      %489 = bitcast {}* %486 to {}**
      store {}* inttoptr (i64 4725386192 to {}*), {}** %489, align 8
      call void @ijl_throw({}* %486)
      unreachable

L1306:                                            ; preds = %idxend35
; │└└└
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:953
; ││┌ @ tuple.jl:398 within `==`
; │││┌ @ tuple.jl:402 within `_eq`
; ││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
       %.not288.not = icmp eq i64 %738, %741
; ││└└└
    br i1 %.not288.not, label %L1322, label %L1325

L1322:                                            ; preds = %L1306
    store {}* %737, {}** %28, align 16
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954 @ array.jl:346 @ array.jl:322
    %490 = call nonnull {}* @"j__copyto_impl!_2830"({}* nonnull %319, i64 signext 1, {}* nonnull %737, i64 signext 1, i64 signext %738) #0
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954
    br label %L1421

L1325:                                            ; preds = %L1306
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:934 within `broadcast_unalias`
        %.not323 = icmp eq {}* %.pre-phi.pre-phi.pre-phi, %737
        br i1 %.not323, label %L1355, label %L1328

L1328:                                            ; preds = %L1325
; ││││││┌ @ abstractarray.jl:1427 within `unalias`
; │││││││┌ @ abstractarray.jl:1462 within `mightalias`
; ││││││││┌ @ reflection.jl:593 within `isbits`
; │││││││││┌ @ Base.jl:33 within `getproperty`
            %491 = load i8, i8* inttoptr (i64 4708472424 to i8*), align 8
; ││││││││└└
          %492 = and i8 %491, 8
          %.not327.not = icmp eq i8 %492, 0
          br i1 %.not327.not, label %L1338, label %L1355

L1338:                                            ; preds = %L1328
; ││││││││┌ @ abstractarray.jl:1486 within `dataids`
; │││││││││┌ @ abstractarray.jl:1187 within `pointer`
; ││││││││││┌ @ pointer.jl:65 within `unsafe_convert`
             %493 = bitcast {}* %319 to i8**
             %494 = load i8*, i8** %493, align 8
             %495 = bitcast {}* %737 to i8**
             %496 = load i8*, i8** %495, align 8
; ││││││││└└└
; ││││││││┌ @ abstractarray.jl:1469 within `_isdisjoint`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %497 = icmp eq i8* %494, %496
; │││││││└└└└
         br i1 %497, label %L1350, label %L1355

L1350:                                            ; preds = %L1338
         store {}* %737, {}** %28, align 16
; │││││││┌ @ abstractarray.jl:1443 within `unaliascopy`
; ││││││││┌ @ array.jl:369 within `copy`
           %498 = call nonnull {}* inttoptr (i64 4320602180 to {}* ({}*)*)({}* nonnull %737)
; └└└└└└└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:39 within `advance_check2`
; ┌ @ reducedim.jl:1007 within `all`
; │┌ @ reducedim.jl:1007 within `#all#795`
; ││┌ @ reduce.jl:1260 within `_all`
; │││┌ @ array.jl:898 within `iterate`
      br label %L1355

L1355:                                            ; preds = %L1350, %L1338, %L1328, %L1325
      %value_phi118 = phi {}* [ %737, %L1325 ], [ %498, %L1350 ], [ %737, %L1338 ], [ %737, %L1328 ]
; └└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:28 within `advance_check2`
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:72 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %.not324 = icmp eq i64 %738, 0
; │││└
     br i1 %.not324, label %L1421, label %L1376.lr.ph

L1376.lr.ph:                                      ; preds = %L1355
; ││└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:637 within `extrude`
; ││││││┌ @ broadcast.jl:586 within `newindexer`
; │││││││┌ @ abstractarray.jl:95 within `axes`
; ││││││││┌ @ array.jl:151 within `size`
           %499 = bitcast {}* %value_phi118 to { i8*, i64, i16, i16, i32 }*
           %500 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %499, i64 0, i32 1
           %501 = load i64, i64* %500, align 8
; │││││││└└
; │││││││┌ @ broadcast.jl:587 within `shapeindexer`
; ││││││││┌ @ broadcast.jl:592 within `_newindexer`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %.not326 = icmp eq i64 %501, 1
             %502 = bitcast {}* %value_phi118 to double**
             %503 = load double*, double** %502, align 8
             %504 = load double*, double** %384, align 8
; ││└└└└└└└└└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:75 within `macro expansion`
     %min.iters.check1476 = icmp ult i64 %738, 4
     br i1 %.not326, label %L1376.us.preheader, label %L1376.preheader

L1376.preheader:                                  ; preds = %L1376.lr.ph
     br i1 %min.iters.check1476, label %L1376, label %vector.memcheck1488

vector.memcheck1488:                              ; preds = %L1376.preheader
     %scevgep1489 = getelementptr double, double* %504, i64 %738
     %scevgep1491 = getelementptr double, double* %503, i64 %738
     %bound01493 = icmp ult double* %504, %scevgep1491
     %bound11494 = icmp ult double* %503, %scevgep1489
     %found.conflict1495 = and i1 %bound01493, %bound11494
     br i1 %found.conflict1495, label %L1376, label %vector.ph1501

vector.ph1501:                                    ; preds = %vector.memcheck1488
     %n.vec1503 = and i64 %738, 9223372036854775804
     br label %vector.body1499

vector.body1499:                                  ; preds = %vector.body1499, %vector.ph1501
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1504 = phi i64 [ 0, %vector.ph1501 ], [ %index.next1505, %vector.body1499 ]
      %505 = getelementptr inbounds double, double* %503, i64 %index1504
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %506 = bitcast double* %505 to <2 x double>*
          %wide.load1508 = load <2 x double>, <2 x double>* %506, align 8
          %507 = getelementptr inbounds double, double* %505, i64 2
          %508 = bitcast double* %507 to <2 x double>*
          %wide.load1509 = load <2 x double>, <2 x double>* %508, align 8
; │││└└└└└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %509 = getelementptr inbounds double, double* %504, i64 %index1504
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ array.jl:966 within `setindex!`
      %510 = bitcast double* %509 to <2 x double>*
      store <2 x double> %wide.load1508, <2 x double>* %510, align 8
      %511 = getelementptr inbounds double, double* %509, i64 2
      %512 = bitcast double* %511 to <2 x double>*
      store <2 x double> %wide.load1509, <2 x double>* %512, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index.next1505 = add nuw i64 %index1504, 4
      %513 = icmp eq i64 %index.next1505, %n.vec1503
      br i1 %513, label %middle.block1497, label %vector.body1499

middle.block1497:                                 ; preds = %vector.body1499
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1507 = icmp eq i64 %738, %n.vec1503
     br i1 %cmp.n1507, label %L1421, label %L1376

L1376.us.preheader:                               ; preds = %L1376.lr.ph
     br i1 %min.iters.check1476, label %L1376.us, label %vector.memcheck1464

vector.memcheck1464:                              ; preds = %L1376.us.preheader
     %scevgep1465 = getelementptr double, double* %504, i64 %738
     %scevgep1467 = getelementptr double, double* %503, i64 1
     %bound01469 = icmp ult double* %504, %scevgep1467
     %bound11470 = icmp ult double* %503, %scevgep1465
     %found.conflict1471 = and i1 %bound01469, %bound11470
     br i1 %found.conflict1471, label %L1376.us, label %vector.ph1477

vector.ph1477:                                    ; preds = %vector.memcheck1464
     %n.vec1479 = and i64 %738, 9223372036854775804
     br label %vector.body1475

vector.body1475:                                  ; preds = %vector.body1475, %vector.ph1477
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1480 = phi i64 [ 0, %vector.ph1477 ], [ %index.next1481, %vector.body1475 ]
      %514 = load double, double* %503, align 8
      %broadcast.splatinsert1484 = insertelement <2 x double> poison, double %514, i32 0
      %broadcast.splat1485 = shufflevector <2 x double> %broadcast.splatinsert1484, <2 x double> poison, <2 x i32> zeroinitializer
      %515 = getelementptr inbounds double, double* %504, i64 %index1480
      %516 = bitcast double* %515 to <2 x double>*
      store <2 x double> %broadcast.splat1485, <2 x double>* %516, align 8
      %517 = getelementptr inbounds double, double* %515, i64 2
      %518 = bitcast double* %517 to <2 x double>*
      store <2 x double> %broadcast.splat1485, <2 x double>* %518, align 8
      %index.next1481 = add nuw i64 %index1480, 4
      %519 = icmp eq i64 %index.next1481, %n.vec1479
      br i1 %519, label %middle.block1473, label %vector.body1475

middle.block1473:                                 ; preds = %vector.body1475
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1483 = icmp eq i64 %738, %n.vec1479
     br i1 %cmp.n1483, label %L1421, label %L1376.us

L1376.us:                                         ; preds = %L1376.us, %middle.block1473, %vector.memcheck1464, %L1376.us.preheader
     %value_phi119461.us = phi i64 [ %522, %L1376.us ], [ %n.vec1479, %middle.block1473 ], [ 0, %L1376.us.preheader ], [ 0, %vector.memcheck1464 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %520 = load double, double* %503, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %521 = getelementptr inbounds double, double* %504, i64 %value_phi119461.us
      store double %520, double* %521, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %522 = add nuw nsw i64 %value_phi119461.us, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond834.not = icmp eq i64 %522, %738
; │││└
     br i1 %exitcond834.not, label %L1421, label %L1376.us

L1376:                                            ; preds = %L1376, %middle.block1497, %vector.memcheck1488, %L1376.preheader
     %value_phi119461 = phi i64 [ %526, %L1376 ], [ %n.vec1503, %middle.block1497 ], [ 0, %L1376.preheader ], [ 0, %vector.memcheck1488 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %523 = getelementptr inbounds double, double* %503, i64 %value_phi119461
          %524 = load double, double* %523, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %525 = getelementptr inbounds double, double* %504, i64 %value_phi119461
      store double %524, double* %525, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %526 = add nuw nsw i64 %value_phi119461, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond828.not = icmp eq i64 %526, %738
; │││└
     br i1 %exitcond828.not, label %L1421, label %L1376

L1421:                                            ; preds = %L1376, %L1376.us, %middle.block1473, %middle.block1497, %L1355, %L1322, %L944
; └└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:31 within `advance_check2`
; ┌ @ array.jl:924 within `getindex`
   %527 = load i64, i64* %386, align 8
   %528 = icmp ult i64 %.pre850, %527
   br i1 %528, label %idxend40, label %oob39

L1446:                                            ; preds = %idxend40
; └
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ broadcast.jl:284 within `instantiate`
; ││┌ @ broadcast.jl:543 within `check_broadcast_axes`
; │││┌ @ broadcast.jl:540 within `check_broadcast_shape`
      %ptls_field16191632 = getelementptr inbounds {}**, {}*** %9, i64 2
      %529 = bitcast {}*** %ptls_field16191632 to i8**
      %ptls_load162016331634 = load i8*, i8** %529, align 8
      %530 = call noalias nonnull {}* @ijl_gc_pool_alloc(i8* %ptls_load162016331634, i32 1392, i32 16) #7
      %531 = bitcast {}* %530 to i64*
      %532 = getelementptr inbounds i64, i64* %531, i64 -1
      store atomic i64 4712015424, i64* %532 unordered, align 8
      %533 = bitcast {}* %530 to {}**
      store {}* inttoptr (i64 4725386192 to {}*), {}** %533, align 8
      call void @ijl_throw({}* %530)
      unreachable

L1465:                                            ; preds = %idxend40
; │└└└
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:953
; ││┌ @ tuple.jl:398 within `==`
; │││┌ @ tuple.jl:402 within `_eq`
; ││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
       %.not290.not = icmp eq i64 %752, %755
; ││└└└
    br i1 %.not290.not, label %L1481, label %L1484

L1481:                                            ; preds = %L1465
    store {}* %751, {}** %28, align 16
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954 @ array.jl:346 @ array.jl:322
    %534 = call nonnull {}* @"j__copyto_impl!_2833"({}* nonnull %43, i64 signext 1, {}* nonnull %751, i64 signext 1, i64 signext %752) #0
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954
    br label %L1579

L1484:                                            ; preds = %L1465
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:934 within `broadcast_unalias`
        %.not316 = icmp eq {}* %43, %751
        br i1 %.not316, label %L1514, label %L1487

L1487:                                            ; preds = %L1484
; ││││││┌ @ abstractarray.jl:1427 within `unalias`
; │││││││┌ @ abstractarray.jl:1462 within `mightalias`
; ││││││││┌ @ reflection.jl:593 within `isbits`
; │││││││││┌ @ Base.jl:33 within `getproperty`
            %535 = load i8, i8* inttoptr (i64 4708472424 to i8*), align 8
; ││││││││└└
          %536 = and i8 %535, 8
          %.not320.not = icmp eq i8 %536, 0
          br i1 %.not320.not, label %L1497, label %L1514

L1497:                                            ; preds = %L1487
; ││││││││┌ @ abstractarray.jl:1486 within `dataids`
; │││││││││┌ @ abstractarray.jl:1187 within `pointer`
; ││││││││││┌ @ pointer.jl:65 within `unsafe_convert`
             %537 = bitcast {}* %43 to i8**
             %538 = load i8*, i8** %537, align 8
             %539 = bitcast {}* %751 to i8**
             %540 = load i8*, i8** %539, align 8
; ││││││││└└└
; ││││││││┌ @ abstractarray.jl:1469 within `_isdisjoint`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %541 = icmp eq i8* %538, %540
; │││││││└└└└
         br i1 %541, label %L1509, label %L1514

L1509:                                            ; preds = %L1497
         store {}* %751, {}** %28, align 16
; │││││││┌ @ abstractarray.jl:1443 within `unaliascopy`
; ││││││││┌ @ array.jl:369 within `copy`
           %542 = call nonnull {}* inttoptr (i64 4320602180 to {}* ({}*)*)({}* nonnull %751)
; └└└└└└└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:39 within `advance_check2`
; ┌ @ reducedim.jl:1007 within `all`
; │┌ @ reducedim.jl:1007 within `#all#795`
; ││┌ @ reduce.jl:1260 within `_all`
; │││┌ @ array.jl:898 within `iterate`
      br label %L1514

L1514:                                            ; preds = %L1509, %L1497, %L1487, %L1484
      %value_phi108 = phi {}* [ %751, %L1484 ], [ %542, %L1509 ], [ %751, %L1497 ], [ %751, %L1487 ]
; └└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:31 within `advance_check2`
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:72 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %.not317 = icmp eq i64 %752, 0
; │││└
     br i1 %.not317, label %L1579, label %L1535.lr.ph

L1535.lr.ph:                                      ; preds = %L1514
; ││└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:637 within `extrude`
; ││││││┌ @ broadcast.jl:586 within `newindexer`
; │││││││┌ @ abstractarray.jl:95 within `axes`
; ││││││││┌ @ array.jl:151 within `size`
           %543 = bitcast {}* %value_phi108 to { i8*, i64, i16, i16, i32 }*
           %544 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %543, i64 0, i32 1
           %545 = load i64, i64* %544, align 8
; │││││││└└
; │││││││┌ @ broadcast.jl:587 within `shapeindexer`
; ││││││││┌ @ broadcast.jl:592 within `_newindexer`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %.not319 = icmp eq i64 %545, 1
             %546 = bitcast {}* %value_phi108 to double**
             %547 = load double*, double** %546, align 8
             %548 = load double*, double** %387, align 8
; ││└└└└└└└└└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:75 within `macro expansion`
     %min.iters.check1430 = icmp ult i64 %752, 4
     br i1 %.not319, label %L1535.us.preheader, label %L1535.preheader

L1535.preheader:                                  ; preds = %L1535.lr.ph
     br i1 %min.iters.check1430, label %L1535, label %vector.memcheck1442

vector.memcheck1442:                              ; preds = %L1535.preheader
     %scevgep1443 = getelementptr double, double* %548, i64 %752
     %scevgep1445 = getelementptr double, double* %547, i64 %752
     %bound01447 = icmp ult double* %548, %scevgep1445
     %bound11448 = icmp ult double* %547, %scevgep1443
     %found.conflict1449 = and i1 %bound01447, %bound11448
     br i1 %found.conflict1449, label %L1535, label %vector.ph1455

vector.ph1455:                                    ; preds = %vector.memcheck1442
     %n.vec1457 = and i64 %752, 9223372036854775804
     br label %vector.body1453

vector.body1453:                                  ; preds = %vector.body1453, %vector.ph1455
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1458 = phi i64 [ 0, %vector.ph1455 ], [ %index.next1459, %vector.body1453 ]
      %549 = getelementptr inbounds double, double* %547, i64 %index1458
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %550 = bitcast double* %549 to <2 x double>*
          %wide.load1462 = load <2 x double>, <2 x double>* %550, align 8
          %551 = getelementptr inbounds double, double* %549, i64 2
          %552 = bitcast double* %551 to <2 x double>*
          %wide.load1463 = load <2 x double>, <2 x double>* %552, align 8
; │││└└└└└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %553 = getelementptr inbounds double, double* %548, i64 %index1458
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ array.jl:966 within `setindex!`
      %554 = bitcast double* %553 to <2 x double>*
      store <2 x double> %wide.load1462, <2 x double>* %554, align 8
      %555 = getelementptr inbounds double, double* %553, i64 2
      %556 = bitcast double* %555 to <2 x double>*
      store <2 x double> %wide.load1463, <2 x double>* %556, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index.next1459 = add nuw i64 %index1458, 4
      %557 = icmp eq i64 %index.next1459, %n.vec1457
      br i1 %557, label %middle.block1451, label %vector.body1453

middle.block1451:                                 ; preds = %vector.body1453
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1461 = icmp eq i64 %752, %n.vec1457
     br i1 %cmp.n1461, label %L1579, label %L1535

L1535.us.preheader:                               ; preds = %L1535.lr.ph
     br i1 %min.iters.check1430, label %L1535.us, label %vector.memcheck1418

vector.memcheck1418:                              ; preds = %L1535.us.preheader
     %scevgep1419 = getelementptr double, double* %548, i64 %752
     %scevgep1421 = getelementptr double, double* %547, i64 1
     %bound01423 = icmp ult double* %548, %scevgep1421
     %bound11424 = icmp ult double* %547, %scevgep1419
     %found.conflict1425 = and i1 %bound01423, %bound11424
     br i1 %found.conflict1425, label %L1535.us, label %vector.ph1431

vector.ph1431:                                    ; preds = %vector.memcheck1418
     %n.vec1433 = and i64 %752, 9223372036854775804
     br label %vector.body1429

vector.body1429:                                  ; preds = %vector.body1429, %vector.ph1431
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1434 = phi i64 [ 0, %vector.ph1431 ], [ %index.next1435, %vector.body1429 ]
      %558 = load double, double* %547, align 8
      %broadcast.splatinsert1438 = insertelement <2 x double> poison, double %558, i32 0
      %broadcast.splat1439 = shufflevector <2 x double> %broadcast.splatinsert1438, <2 x double> poison, <2 x i32> zeroinitializer
      %559 = getelementptr inbounds double, double* %548, i64 %index1434
      %560 = bitcast double* %559 to <2 x double>*
      store <2 x double> %broadcast.splat1439, <2 x double>* %560, align 8
      %561 = getelementptr inbounds double, double* %559, i64 2
      %562 = bitcast double* %561 to <2 x double>*
      store <2 x double> %broadcast.splat1439, <2 x double>* %562, align 8
      %index.next1435 = add nuw i64 %index1434, 4
      %563 = icmp eq i64 %index.next1435, %n.vec1433
      br i1 %563, label %middle.block1427, label %vector.body1429

middle.block1427:                                 ; preds = %vector.body1429
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1437 = icmp eq i64 %752, %n.vec1433
     br i1 %cmp.n1437, label %L1579, label %L1535.us

L1535.us:                                         ; preds = %L1535.us, %middle.block1427, %vector.memcheck1418, %L1535.us.preheader
     %value_phi109463.us = phi i64 [ %566, %L1535.us ], [ %n.vec1433, %middle.block1427 ], [ 0, %L1535.us.preheader ], [ 0, %vector.memcheck1418 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %564 = load double, double* %547, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %565 = getelementptr inbounds double, double* %548, i64 %value_phi109463.us
      store double %564, double* %565, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %566 = add nuw nsw i64 %value_phi109463.us, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond835.not = icmp eq i64 %566, %752
; │││└
     br i1 %exitcond835.not, label %L1579, label %L1535.us

L1535:                                            ; preds = %L1535, %middle.block1451, %vector.memcheck1442, %L1535.preheader
     %value_phi109463 = phi i64 [ %570, %L1535 ], [ %n.vec1457, %middle.block1451 ], [ 0, %L1535.preheader ], [ 0, %vector.memcheck1442 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %567 = getelementptr inbounds double, double* %547, i64 %value_phi109463
          %568 = load double, double* %567, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %569 = getelementptr inbounds double, double* %548, i64 %value_phi109463
      store double %568, double* %569, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %570 = add nuw nsw i64 %value_phi109463, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond829.not = icmp eq i64 %570, %752
; │││└
     br i1 %exitcond829.not, label %L1579, label %L1535

L1579:                                            ; preds = %L1535, %L1535.us, %middle.block1427, %middle.block1451, %L1514, %L1481
; └└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:32 within `advance_check2`
; ┌ @ array.jl:924 within `getindex`
   %571 = load i64, i64* %386, align 8
   %572 = icmp ult i64 %.pre850, %571
   br i1 %572, label %idxend45, label %oob44

L1604:                                            ; preds = %idxend45
; └
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ broadcast.jl:284 within `instantiate`
; ││┌ @ broadcast.jl:543 within `check_broadcast_axes`
; │││┌ @ broadcast.jl:540 within `check_broadcast_shape`
      %ptls_field16211629 = getelementptr inbounds {}**, {}*** %9, i64 2
      %573 = bitcast {}*** %ptls_field16211629 to i8**
      %ptls_load162216301631 = load i8*, i8** %573, align 8
      %574 = call noalias nonnull {}* @ijl_gc_pool_alloc(i8* %ptls_load162216301631, i32 1392, i32 16) #7
      %575 = bitcast {}* %574 to i64*
      %576 = getelementptr inbounds i64, i64* %575, i64 -1
      store atomic i64 4712015424, i64* %576 unordered, align 8
      %577 = bitcast {}* %574 to {}**
      store {}* inttoptr (i64 4725386192 to {}*), {}** %577, align 8
      call void @ijl_throw({}* %574)
      unreachable

L1623:                                            ; preds = %idxend45
; │└└└
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:953
; ││┌ @ tuple.jl:398 within `==`
; │││┌ @ tuple.jl:402 within `_eq`
; ││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
       %.not292.not = icmp eq i64 %766, %769
; ││└└└
    br i1 %.not292.not, label %L1639, label %L1642

L1639:                                            ; preds = %L1623
    store {}* %765, {}** %28, align 16
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954 @ array.jl:346 @ array.jl:322
    %578 = call nonnull {}* @"j__copyto_impl!_2836"({}* nonnull %99, i64 signext 1, {}* nonnull %765, i64 signext 1, i64 signext %766) #0
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954
    br label %L1737

L1642:                                            ; preds = %L1623
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:934 within `broadcast_unalias`
        %.not309 = icmp eq {}* %99, %765
        br i1 %.not309, label %L1672, label %L1645

L1645:                                            ; preds = %L1642
; ││││││┌ @ abstractarray.jl:1427 within `unalias`
; │││││││┌ @ abstractarray.jl:1462 within `mightalias`
; ││││││││┌ @ reflection.jl:593 within `isbits`
; │││││││││┌ @ Base.jl:33 within `getproperty`
            %579 = load i8, i8* inttoptr (i64 4708472424 to i8*), align 8
; ││││││││└└
          %580 = and i8 %579, 8
          %.not313.not = icmp eq i8 %580, 0
          br i1 %.not313.not, label %L1655, label %L1672

L1655:                                            ; preds = %L1645
; ││││││││┌ @ abstractarray.jl:1486 within `dataids`
; │││││││││┌ @ abstractarray.jl:1187 within `pointer`
; ││││││││││┌ @ pointer.jl:65 within `unsafe_convert`
             %581 = bitcast {}* %99 to i8**
             %582 = load i8*, i8** %581, align 8
             %583 = bitcast {}* %765 to i8**
             %584 = load i8*, i8** %583, align 8
; ││││││││└└└
; ││││││││┌ @ abstractarray.jl:1469 within `_isdisjoint`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %585 = icmp eq i8* %582, %584
; │││││││└└└└
         br i1 %585, label %L1667, label %L1672

L1667:                                            ; preds = %L1655
         store {}* %765, {}** %28, align 16
; │││││││┌ @ abstractarray.jl:1443 within `unaliascopy`
; ││││││││┌ @ array.jl:369 within `copy`
           %586 = call nonnull {}* inttoptr (i64 4320602180 to {}* ({}*)*)({}* nonnull %765)
; └└└└└└└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:39 within `advance_check2`
; ┌ @ reducedim.jl:1007 within `all`
; │┌ @ reducedim.jl:1007 within `#all#795`
; ││┌ @ reduce.jl:1260 within `_all`
; │││┌ @ array.jl:898 within `iterate`
      br label %L1672

L1672:                                            ; preds = %L1667, %L1655, %L1645, %L1642
      %value_phi98 = phi {}* [ %765, %L1642 ], [ %586, %L1667 ], [ %765, %L1655 ], [ %765, %L1645 ]
; └└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:32 within `advance_check2`
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:72 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %.not310 = icmp eq i64 %766, 0
; │││└
     br i1 %.not310, label %L1737, label %L1693.lr.ph

L1693.lr.ph:                                      ; preds = %L1672
; ││└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:637 within `extrude`
; ││││││┌ @ broadcast.jl:586 within `newindexer`
; │││││││┌ @ abstractarray.jl:95 within `axes`
; ││││││││┌ @ array.jl:151 within `size`
           %587 = bitcast {}* %value_phi98 to { i8*, i64, i16, i16, i32 }*
           %588 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %587, i64 0, i32 1
           %589 = load i64, i64* %588, align 8
; │││││││└└
; │││││││┌ @ broadcast.jl:587 within `shapeindexer`
; ││││││││┌ @ broadcast.jl:592 within `_newindexer`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %.not312 = icmp eq i64 %589, 1
             %590 = bitcast {}* %value_phi98 to double**
             %591 = load double*, double** %590, align 8
             %592 = load double*, double** %388, align 8
; ││└└└└└└└└└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:75 within `macro expansion`
     %min.iters.check1384 = icmp ult i64 %766, 4
     br i1 %.not312, label %L1693.us.preheader, label %L1693.preheader

L1693.preheader:                                  ; preds = %L1693.lr.ph
     br i1 %min.iters.check1384, label %L1693, label %vector.memcheck1396

vector.memcheck1396:                              ; preds = %L1693.preheader
     %scevgep1397 = getelementptr double, double* %592, i64 %766
     %scevgep1399 = getelementptr double, double* %591, i64 %766
     %bound01401 = icmp ult double* %592, %scevgep1399
     %bound11402 = icmp ult double* %591, %scevgep1397
     %found.conflict1403 = and i1 %bound01401, %bound11402
     br i1 %found.conflict1403, label %L1693, label %vector.ph1409

vector.ph1409:                                    ; preds = %vector.memcheck1396
     %n.vec1411 = and i64 %766, 9223372036854775804
     br label %vector.body1407

vector.body1407:                                  ; preds = %vector.body1407, %vector.ph1409
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1412 = phi i64 [ 0, %vector.ph1409 ], [ %index.next1413, %vector.body1407 ]
      %593 = getelementptr inbounds double, double* %591, i64 %index1412
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %594 = bitcast double* %593 to <2 x double>*
          %wide.load1416 = load <2 x double>, <2 x double>* %594, align 8
          %595 = getelementptr inbounds double, double* %593, i64 2
          %596 = bitcast double* %595 to <2 x double>*
          %wide.load1417 = load <2 x double>, <2 x double>* %596, align 8
; │││└└└└└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %597 = getelementptr inbounds double, double* %592, i64 %index1412
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ array.jl:966 within `setindex!`
      %598 = bitcast double* %597 to <2 x double>*
      store <2 x double> %wide.load1416, <2 x double>* %598, align 8
      %599 = getelementptr inbounds double, double* %597, i64 2
      %600 = bitcast double* %599 to <2 x double>*
      store <2 x double> %wide.load1417, <2 x double>* %600, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index.next1413 = add nuw i64 %index1412, 4
      %601 = icmp eq i64 %index.next1413, %n.vec1411
      br i1 %601, label %middle.block1405, label %vector.body1407

middle.block1405:                                 ; preds = %vector.body1407
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1415 = icmp eq i64 %766, %n.vec1411
     br i1 %cmp.n1415, label %L1737, label %L1693

L1693.us.preheader:                               ; preds = %L1693.lr.ph
     br i1 %min.iters.check1384, label %L1693.us, label %vector.memcheck1372

vector.memcheck1372:                              ; preds = %L1693.us.preheader
     %scevgep1373 = getelementptr double, double* %592, i64 %766
     %scevgep1375 = getelementptr double, double* %591, i64 1
     %bound01377 = icmp ult double* %592, %scevgep1375
     %bound11378 = icmp ult double* %591, %scevgep1373
     %found.conflict1379 = and i1 %bound01377, %bound11378
     br i1 %found.conflict1379, label %L1693.us, label %vector.ph1385

vector.ph1385:                                    ; preds = %vector.memcheck1372
     %n.vec1387 = and i64 %766, 9223372036854775804
     br label %vector.body1383

vector.body1383:                                  ; preds = %vector.body1383, %vector.ph1385
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1388 = phi i64 [ 0, %vector.ph1385 ], [ %index.next1389, %vector.body1383 ]
      %602 = load double, double* %591, align 8
      %broadcast.splatinsert1392 = insertelement <2 x double> poison, double %602, i32 0
      %broadcast.splat1393 = shufflevector <2 x double> %broadcast.splatinsert1392, <2 x double> poison, <2 x i32> zeroinitializer
      %603 = getelementptr inbounds double, double* %592, i64 %index1388
      %604 = bitcast double* %603 to <2 x double>*
      store <2 x double> %broadcast.splat1393, <2 x double>* %604, align 8
      %605 = getelementptr inbounds double, double* %603, i64 2
      %606 = bitcast double* %605 to <2 x double>*
      store <2 x double> %broadcast.splat1393, <2 x double>* %606, align 8
      %index.next1389 = add nuw i64 %index1388, 4
      %607 = icmp eq i64 %index.next1389, %n.vec1387
      br i1 %607, label %middle.block1381, label %vector.body1383

middle.block1381:                                 ; preds = %vector.body1383
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1391 = icmp eq i64 %766, %n.vec1387
     br i1 %cmp.n1391, label %L1737, label %L1693.us

L1693.us:                                         ; preds = %L1693.us, %middle.block1381, %vector.memcheck1372, %L1693.us.preheader
     %value_phi99465.us = phi i64 [ %610, %L1693.us ], [ %n.vec1387, %middle.block1381 ], [ 0, %L1693.us.preheader ], [ 0, %vector.memcheck1372 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %608 = load double, double* %591, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %609 = getelementptr inbounds double, double* %592, i64 %value_phi99465.us
      store double %608, double* %609, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %610 = add nuw nsw i64 %value_phi99465.us, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond836.not = icmp eq i64 %610, %766
; │││└
     br i1 %exitcond836.not, label %L1737, label %L1693.us

L1693:                                            ; preds = %L1693, %middle.block1405, %vector.memcheck1396, %L1693.preheader
     %value_phi99465 = phi i64 [ %614, %L1693 ], [ %n.vec1411, %middle.block1405 ], [ 0, %L1693.preheader ], [ 0, %vector.memcheck1396 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %611 = getelementptr inbounds double, double* %591, i64 %value_phi99465
          %612 = load double, double* %611, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %613 = getelementptr inbounds double, double* %592, i64 %value_phi99465
      store double %612, double* %613, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %614 = add nuw nsw i64 %value_phi99465, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond830.not = icmp eq i64 %614, %766
; │││└
     br i1 %exitcond830.not, label %L1737, label %L1693

L1737:                                            ; preds = %L1693, %L1693.us, %middle.block1381, %middle.block1405, %L1672, %L1639
; └└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:33 within `advance_check2`
; ┌ @ array.jl:924 within `getindex`
   %615 = load i64, i64* %386, align 8
   %616 = icmp ult i64 %.pre850, %615
   br i1 %616, label %idxend50, label %oob49

L1762:                                            ; preds = %idxend50
; └
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ broadcast.jl:284 within `instantiate`
; ││┌ @ broadcast.jl:543 within `check_broadcast_axes`
; │││┌ @ broadcast.jl:540 within `check_broadcast_shape`
      %ptls_field16231626 = getelementptr inbounds {}**, {}*** %9, i64 2
      %617 = bitcast {}*** %ptls_field16231626 to i8**
      %ptls_load162416271628 = load i8*, i8** %617, align 8
      %618 = call noalias nonnull {}* @ijl_gc_pool_alloc(i8* %ptls_load162416271628, i32 1392, i32 16) #7
      %619 = bitcast {}* %618 to i64*
      %620 = getelementptr inbounds i64, i64* %619, i64 -1
      store atomic i64 4712015424, i64* %620 unordered, align 8
      %621 = bitcast {}* %618 to {}**
      store {}* inttoptr (i64 4725386192 to {}*), {}** %621, align 8
      call void @ijl_throw({}* %618)
      unreachable

L1781:                                            ; preds = %idxend50
; │└└└
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:953
; ││┌ @ tuple.jl:398 within `==`
; │││┌ @ tuple.jl:402 within `_eq`
; ││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
       %.not294.not = icmp eq i64 %780, %783
; ││└└└
    br i1 %.not294.not, label %L1797, label %L1800

L1797:                                            ; preds = %L1781
    store {}* %779, {}** %28, align 16
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954 @ array.jl:346 @ array.jl:322
    %622 = call nonnull {}* @"j__copyto_impl!_2839"({}* nonnull %155, i64 signext 1, {}* nonnull %779, i64 signext 1, i64 signext %780) #0
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:954
    br label %L1895

L1800:                                            ; preds = %L1781
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:934 within `broadcast_unalias`
        %.not302 = icmp eq {}* %155, %779
        br i1 %.not302, label %L1830, label %L1803

L1803:                                            ; preds = %L1800
; ││││││┌ @ abstractarray.jl:1427 within `unalias`
; │││││││┌ @ abstractarray.jl:1462 within `mightalias`
; ││││││││┌ @ reflection.jl:593 within `isbits`
; │││││││││┌ @ Base.jl:33 within `getproperty`
            %623 = load i8, i8* inttoptr (i64 4708472424 to i8*), align 8
; ││││││││└└
          %624 = and i8 %623, 8
          %.not306.not = icmp eq i8 %624, 0
          br i1 %.not306.not, label %L1813, label %L1830

L1813:                                            ; preds = %L1803
; ││││││││┌ @ abstractarray.jl:1486 within `dataids`
; │││││││││┌ @ abstractarray.jl:1187 within `pointer`
; ││││││││││┌ @ pointer.jl:65 within `unsafe_convert`
             %625 = bitcast {}* %155 to i8**
             %626 = load i8*, i8** %625, align 8
             %627 = bitcast {}* %779 to i8**
             %628 = load i8*, i8** %627, align 8
; ││││││││└└└
; ││││││││┌ @ abstractarray.jl:1469 within `_isdisjoint`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %629 = icmp eq i8* %626, %628
; │││││││└└└└
         br i1 %629, label %L1825, label %L1830

L1825:                                            ; preds = %L1813
         store {}* %779, {}** %28, align 16
; │││││││┌ @ abstractarray.jl:1443 within `unaliascopy`
; ││││││││┌ @ array.jl:369 within `copy`
           %630 = call nonnull {}* inttoptr (i64 4320602180 to {}* ({}*)*)({}* nonnull %779)
; └└└└└└└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:39 within `advance_check2`
; ┌ @ reducedim.jl:1007 within `all`
; │┌ @ reducedim.jl:1007 within `#all#795`
; ││┌ @ reduce.jl:1260 within `_all`
; │││┌ @ array.jl:898 within `iterate`
      br label %L1830

L1830:                                            ; preds = %L1825, %L1813, %L1803, %L1800
      %value_phi88 = phi {}* [ %779, %L1800 ], [ %630, %L1825 ], [ %779, %L1813 ], [ %779, %L1803 ]
; └└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:33 within `advance_check2`
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:72 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %.not303 = icmp eq i64 %780, 0
; │││└
     br i1 %.not303, label %L1895, label %L1851.lr.ph

L1851.lr.ph:                                      ; preds = %L1830
; ││└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:957
; ││┌ @ broadcast.jl:940 within `preprocess`
; │││┌ @ broadcast.jl:944 within `preprocess_args`
; ││││┌ @ broadcast.jl:941 within `preprocess`
; │││││┌ @ broadcast.jl:637 within `extrude`
; ││││││┌ @ broadcast.jl:586 within `newindexer`
; │││││││┌ @ abstractarray.jl:95 within `axes`
; ││││││││┌ @ array.jl:151 within `size`
           %631 = bitcast {}* %value_phi88 to { i8*, i64, i16, i16, i32 }*
           %632 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %631, i64 0, i32 1
           %633 = load i64, i64* %632, align 8
; │││││││└└
; │││││││┌ @ broadcast.jl:587 within `shapeindexer`
; ││││││││┌ @ broadcast.jl:592 within `_newindexer`
; │││││││││┌ @ operators.jl:282 within `!=`
; ││││││││││┌ @ promotion.jl:477 within `==`
             %.not305 = icmp eq i64 %633, 1
             %634 = bitcast {}* %value_phi88 to double**
             %635 = load double*, double** %634, align 8
             %636 = load double*, double** %389, align 8
; ││└└└└└└└└└
; ││ @ broadcast.jl:913 within `copyto!` @ broadcast.jl:960
; ││┌ @ simdloop.jl:75 within `macro expansion`
     %min.iters.check1338 = icmp ult i64 %780, 4
     br i1 %.not305, label %L1851.us.preheader, label %L1851.preheader

L1851.preheader:                                  ; preds = %L1851.lr.ph
     br i1 %min.iters.check1338, label %L1851, label %vector.memcheck1350

vector.memcheck1350:                              ; preds = %L1851.preheader
     %scevgep1351 = getelementptr double, double* %636, i64 %780
     %scevgep1353 = getelementptr double, double* %635, i64 %780
     %bound01355 = icmp ult double* %636, %scevgep1353
     %bound11356 = icmp ult double* %635, %scevgep1351
     %found.conflict1357 = and i1 %bound01355, %bound11356
     br i1 %found.conflict1357, label %L1851, label %vector.ph1363

vector.ph1363:                                    ; preds = %vector.memcheck1350
     %n.vec1365 = and i64 %780, 9223372036854775804
     br label %vector.body1361

vector.body1361:                                  ; preds = %vector.body1361, %vector.ph1363
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1366 = phi i64 [ 0, %vector.ph1363 ], [ %index.next1367, %vector.body1361 ]
      %637 = getelementptr inbounds double, double* %635, i64 %index1366
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %638 = bitcast double* %637 to <2 x double>*
          %wide.load1370 = load <2 x double>, <2 x double>* %638, align 8
          %639 = getelementptr inbounds double, double* %637, i64 2
          %640 = bitcast double* %639 to <2 x double>*
          %wide.load1371 = load <2 x double>, <2 x double>* %640, align 8
; │││└└└└└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %641 = getelementptr inbounds double, double* %636, i64 %index1366
; │││└
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ array.jl:966 within `setindex!`
      %642 = bitcast double* %641 to <2 x double>*
      store <2 x double> %wide.load1370, <2 x double>* %642, align 8
      %643 = getelementptr inbounds double, double* %641, i64 2
      %644 = bitcast double* %643 to <2 x double>*
      store <2 x double> %wide.load1371, <2 x double>* %644, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index.next1367 = add nuw i64 %index1366, 4
      %645 = icmp eq i64 %index.next1367, %n.vec1365
      br i1 %645, label %middle.block1359, label %vector.body1361

middle.block1359:                                 ; preds = %vector.body1361
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1369 = icmp eq i64 %780, %n.vec1365
     br i1 %cmp.n1369, label %L1895, label %L1851

L1851.us.preheader:                               ; preds = %L1851.lr.ph
     br i1 %min.iters.check1338, label %L1851.us, label %vector.memcheck1326

vector.memcheck1326:                              ; preds = %L1851.us.preheader
     %scevgep1327 = getelementptr double, double* %636, i64 %780
     %scevgep1329 = getelementptr double, double* %635, i64 1
     %bound01331 = icmp ult double* %636, %scevgep1329
     %bound11332 = icmp ult double* %635, %scevgep1327
     %found.conflict1333 = and i1 %bound01331, %bound11332
     br i1 %found.conflict1333, label %L1851.us, label %vector.ph1339

vector.ph1339:                                    ; preds = %vector.memcheck1326
     %n.vec1341 = and i64 %780, 9223372036854775804
     br label %vector.body1337

vector.body1337:                                  ; preds = %vector.body1337, %vector.ph1339
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %index1342 = phi i64 [ 0, %vector.ph1339 ], [ %index.next1343, %vector.body1337 ]
      %646 = load double, double* %635, align 8
      %broadcast.splatinsert1346 = insertelement <2 x double> poison, double %646, i32 0
      %broadcast.splat1347 = shufflevector <2 x double> %broadcast.splatinsert1346, <2 x double> poison, <2 x i32> zeroinitializer
      %647 = getelementptr inbounds double, double* %636, i64 %index1342
      %648 = bitcast double* %647 to <2 x double>*
      store <2 x double> %broadcast.splat1347, <2 x double>* %648, align 8
      %649 = getelementptr inbounds double, double* %647, i64 2
      %650 = bitcast double* %649 to <2 x double>*
      store <2 x double> %broadcast.splat1347, <2 x double>* %650, align 8
      %index.next1343 = add nuw i64 %index1342, 4
      %651 = icmp eq i64 %index.next1343, %n.vec1341
      br i1 %651, label %middle.block1335, label %vector.body1337

middle.block1335:                                 ; preds = %vector.body1337
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
     %cmp.n1345 = icmp eq i64 %780, %n.vec1341
     br i1 %cmp.n1345, label %L1895, label %L1851.us

L1851.us:                                         ; preds = %L1851.us, %middle.block1335, %vector.memcheck1326, %L1851.us.preheader
     %value_phi89467.us = phi i64 [ %654, %L1851.us ], [ %n.vec1341, %middle.block1335 ], [ 0, %L1851.us.preheader ], [ 0, %vector.memcheck1326 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %652 = load double, double* %635, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %653 = getelementptr inbounds double, double* %636, i64 %value_phi89467.us
      store double %652, double* %653, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %654 = add nuw nsw i64 %value_phi89467.us, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond837.not = icmp eq i64 %654, %780
; │││└
     br i1 %exitcond837.not, label %L1895, label %L1851.us

L1851:                                            ; preds = %L1851, %middle.block1359, %vector.memcheck1350, %L1851.preheader
     %value_phi89467 = phi i64 [ %658, %L1851 ], [ %n.vec1365, %middle.block1359 ], [ 0, %L1851.preheader ], [ 0, %vector.memcheck1350 ]
; │││ @ simdloop.jl:77 within `macro expansion` @ broadcast.jl:961
; │││┌ @ broadcast.jl:597 within `getindex`
; ││││┌ @ broadcast.jl:642 within `_broadcast_getindex`
; │││││┌ @ broadcast.jl:667 within `_getindex`
; ││││││┌ @ broadcast.jl:636 within `_broadcast_getindex`
; │││││││┌ @ array.jl:924 within `getindex`
          %655 = getelementptr inbounds double, double* %635, i64 %value_phi89467
          %656 = load double, double* %655, align 8
; │││└└└└└
; │││┌ @ array.jl:966 within `setindex!`
      %657 = getelementptr inbounds double, double* %636, i64 %value_phi89467
      store double %656, double* %657, align 8
; │││└
; │││ @ simdloop.jl:78 within `macro expansion`
; │││┌ @ int.jl:87 within `+`
      %658 = add nuw nsw i64 %value_phi89467, 1
; │││└
; │││ @ simdloop.jl:75 within `macro expansion`
; │││┌ @ int.jl:83 within `<`
      %exitcond831.not = icmp eq i64 %658, %780
; │││└
     br i1 %exitcond831.not, label %L1895, label %L1851

L1895:                                            ; preds = %L1851, %L1851.us, %middle.block1335, %middle.block1359, %L1830, %L1797
; └└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:35 within `advance_check2`
; ┌ @ range.jl:883 within `iterate`
; │┌ @ promotion.jl:477 within `==`
    %.not295.not = icmp eq i64 %value_phi24, 4
; │└
   %659 = add nuw nsw i64 %value_phi24, 1
; └
  br i1 %.not295.not, label %L1908, label %L944

L1908:                                            ; preds = %L1895
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:37 within `advance_check2`
; ┌ @ reducedim.jl:1007 within `all`
; │┌ @ reducedim.jl:1007 within `#all#795`
; ││┌ @ reduce.jl:1250 within `_all`
; │││┌ @ array.jl:898 within `iterate` @ array.jl:898
; ││││┌ @ array.jl:215 within `length`
       %660 = load i64, i64* %45, align 8
; ││││└
; ││││┌ @ int.jl:487 within `<` @ int.jl:480
       %.not296 = icmp eq i64 %660, 0
; ││││└
      br i1 %.not296, label %L1957, label %L1926

L1926:                                            ; preds = %L1908
; ││││┌ @ array.jl:924 within `getindex`
       %661 = load double*, double** %387, align 8
       %662 = load double, double* %661, align 8
; │││└└
; │││ @ reduce.jl:1251 within `_all`
; │││┌ @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:37 within `#79`
; ││││┌ @ float.jl:412 within `<`
       %663 = fcmp uge double %662, 1.000000e+01
; │││└└
; │││ @ reduce.jl:1255 within `_all`
     br i1 %663, label %L2077, label %L1932.lr.ph

L1932.lr.ph:                                      ; preds = %L1926
     %664 = add nuw nsw i64 %660, 1
     br label %L1932

L1932:                                            ; preds = %L1945, %L1932.lr.ph
     %value_phi61455 = phi i64 [ 2, %L1932.lr.ph ], [ %668, %L1945 ]
; │││ @ reduce.jl:1260 within `_all`
; │││┌ @ array.jl:898 within `iterate`
; ││││┌ @ int.jl:487 within `<` @ int.jl:480
       %exitcond825.not = icmp eq i64 %value_phi61455, %664
; ││││└
      br i1 %exitcond825.not, label %L1957, label %L1945

L1945:                                            ; preds = %L1932
; ││││┌ @ int.jl:991 within `-` @ int.jl:86
       %665 = add nsw i64 %value_phi61455, -1
; ││││└
; ││││┌ @ array.jl:924 within `getindex`
       %666 = getelementptr inbounds double, double* %661, i64 %665
       %667 = load double, double* %666, align 8
; ││││└
; ││││┌ @ int.jl:87 within `+`
       %668 = add nuw i64 %value_phi61455, 1
; │││└└
; │││ @ reduce.jl:1251 within `_all`
; │││┌ @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:37 within `#79`
; ││││┌ @ float.jl:412 within `<`
       %669 = fcmp uge double %667, 1.000000e+01
; │││└└
; │││ @ reduce.jl:1255 within `_all`
     br i1 %669, label %L2077, label %L1932

L1957:                                            ; preds = %L1932, %L1908
; └└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:38 within `advance_check2`
; ┌ @ reducedim.jl:1007 within `all`
; │┌ @ reducedim.jl:1007 within `#all#795`
; ││┌ @ reduce.jl:1250 within `_all`
; │││┌ @ array.jl:898 within `iterate` @ array.jl:898
; ││││┌ @ array.jl:215 within `length`
       %670 = load i64, i64* %101, align 8
; ││││└
; ││││┌ @ int.jl:487 within `<` @ int.jl:480
       %.not298 = icmp eq i64 %670, 0
; ││││└
      br i1 %.not298, label %L2007, label %L1976

L1976:                                            ; preds = %L1957
; ││││┌ @ array.jl:924 within `getindex`
       %671 = load double*, double** %388, align 8
       %672 = load double, double* %671, align 8
; │││└└
; │││ @ reduce.jl:1251 within `_all`
; │││┌ @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:38 within `#80`
; ││││┌ @ float.jl:412 within `<`
       %673 = fcmp uge double %672, 1.000000e+01
; │││└└
; │││ @ reduce.jl:1255 within `_all`
     br i1 %673, label %L2074, label %L1982.lr.ph

L1982.lr.ph:                                      ; preds = %L1976
     %674 = add nuw nsw i64 %670, 1
     br label %L1982

L1982:                                            ; preds = %L1995, %L1982.lr.ph
     %value_phi70454 = phi i64 [ 2, %L1982.lr.ph ], [ %678, %L1995 ]
; │││ @ reduce.jl:1260 within `_all`
; │││┌ @ array.jl:898 within `iterate`
; ││││┌ @ int.jl:487 within `<` @ int.jl:480
       %exitcond823.not = icmp eq i64 %value_phi70454, %674
; ││││└
      br i1 %exitcond823.not, label %L2007, label %L1995

L1995:                                            ; preds = %L1982
; ││││┌ @ int.jl:991 within `-` @ int.jl:86
       %675 = add nsw i64 %value_phi70454, -1
; ││││└
; ││││┌ @ array.jl:924 within `getindex`
       %676 = getelementptr inbounds double, double* %671, i64 %675
       %677 = load double, double* %676, align 8
; ││││└
; ││││┌ @ int.jl:87 within `+`
       %678 = add nuw i64 %value_phi70454, 1
; │││└└
; │││ @ reduce.jl:1251 within `_all`
; │││┌ @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:38 within `#80`
; ││││┌ @ float.jl:412 within `<`
       %679 = fcmp uge double %677, 1.000000e+01
; │││└└
; │││ @ reduce.jl:1255 within `_all`
     br i1 %679, label %L2074, label %L1982

L2007:                                            ; preds = %L1982, %L1957
; └└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:39 within `advance_check2`
; ┌ @ reducedim.jl:1007 within `all`
; │┌ @ reducedim.jl:1007 within `#all#795`
; ││┌ @ reduce.jl:1250 within `_all`
; │││┌ @ array.jl:898 within `iterate` @ array.jl:898
; ││││┌ @ array.jl:215 within `length`
       %680 = load i64, i64* %157, align 8
; ││││└
; ││││┌ @ int.jl:487 within `<` @ int.jl:480
       %.not300 = icmp eq i64 %680, 0
; ││││└
      br i1 %.not300, label %L2057, label %L2026

L2026:                                            ; preds = %L2007
; ││││┌ @ array.jl:924 within `getindex`
       %681 = load double*, double** %389, align 8
       %682 = load double, double* %681, align 8
; │││└└
; │││ @ reduce.jl:1251 within `_all`
; │││┌ @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:39 within `#81`
; ││││┌ @ float.jl:412 within `<`
       %683 = fcmp uge double %682, 1.000000e+01
; │││└└
; │││ @ reduce.jl:1255 within `_all`
     br i1 %683, label %L2071, label %L2032.lr.ph

L2032.lr.ph:                                      ; preds = %L2026
     %684 = add nuw nsw i64 %680, 1
     br label %L2032

L2032:                                            ; preds = %L2045, %L2032.lr.ph
     %value_phi79453 = phi i64 [ 2, %L2032.lr.ph ], [ %688, %L2045 ]
; │││ @ reduce.jl:1260 within `_all`
; │││┌ @ array.jl:898 within `iterate`
; ││││┌ @ int.jl:487 within `<` @ int.jl:480
       %exitcond.not = icmp eq i64 %value_phi79453, %684
; ││││└
      br i1 %exitcond.not, label %L2057, label %L2045

L2045:                                            ; preds = %L2032
; ││││┌ @ int.jl:991 within `-` @ int.jl:86
       %685 = add nsw i64 %value_phi79453, -1
; ││││└
; ││││┌ @ array.jl:924 within `getindex`
       %686 = getelementptr inbounds double, double* %681, i64 %685
       %687 = load double, double* %686, align 8
; ││││└
; ││││┌ @ int.jl:87 within `+`
       %688 = add nuw i64 %value_phi79453, 1
; │││└└
; │││ @ reduce.jl:1251 within `_all`
; │││┌ @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:39 within `#81`
; ││││┌ @ float.jl:412 within `<`
       %689 = fcmp uge double %687, 1.000000e+01
; │││└└
; │││ @ reduce.jl:1255 within `_all`
     br i1 %689, label %L2071, label %L2032

L2057:                                            ; preds = %L2032, %L2007
; └└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:41 within `advance_check2`
; ┌ @ Base.jl:38 within `getproperty`
   %690 = load atomic {}*, {}** %26 unordered, align 8
   store {}* %690, {}** %28, align 16
; └
; ┌ @ array.jl:346 within `copyto!` @ array.jl:322
   %691 = call nonnull {}* @"j__copyto_impl!_2840"({}* nonnull %690, i64 signext 1, {}* nonnull %43, i64 signext 1, i64 signext %660) #0
; └
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:42 within `advance_check2`
; ┌ @ Base.jl:38 within `getproperty`
   %692 = load atomic {}*, {}** %33 unordered, align 8
; └
; ┌ @ array.jl:346 within `copyto!`
; │┌ @ array.jl:215 within `length`
    %693 = load i64, i64* %101, align 8
    store {}* %692, {}** %28, align 16
; │└
; │ @ array.jl:346 within `copyto!` @ array.jl:322
   %694 = call nonnull {}* @"j__copyto_impl!_2841"({}* nonnull %692, i64 signext 1, {}* nonnull %99, i64 signext 1, i64 signext %693) #0
; └
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:43 within `advance_check2`
; ┌ @ Base.jl:38 within `getproperty`
   %695 = load atomic {}*, {}** %38 unordered, align 8
; └
; ┌ @ array.jl:346 within `copyto!`
; │┌ @ array.jl:215 within `length`
    %696 = load i64, i64* %157, align 8
    store {}* %695, {}** %28, align 16
; │└
; │ @ array.jl:346 within `copyto!` @ array.jl:322
   %697 = call nonnull {}* @"j__copyto_impl!_2842"({}* nonnull %695, i64 signext 1, {}* nonnull %155, i64 signext 1, i64 signext %696) #0
   %698 = load {}*, {}** %12, align 8
   %699 = bitcast {}*** %9 to {}**
   store {}* %698, {}** %699, align 8
; └
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:45 within `advance_check2`
  ret void

L2071:                                            ; preds = %L2045, %L2026
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:39 within `advance_check2`
  store {}* inttoptr (i64 4446597712 to {}*), {}** %.sub, align 8
  %700 = call nonnull {}* @ijl_apply_generic({}* inttoptr (i64 4710351008 to {}*), {}** nonnull %.sub, i32 1)
  call void @ijl_throw({}* %700)
  unreachable

L2074:                                            ; preds = %L1995, %L1976
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:38 within `advance_check2`
  store {}* inttoptr (i64 4446597328 to {}*), {}** %.sub, align 8
  %701 = call nonnull {}* @ijl_apply_generic({}* inttoptr (i64 4710351008 to {}*), {}** nonnull %.sub, i32 1)
  call void @ijl_throw({}* %701)
  unreachable

L2077:                                            ; preds = %L1945, %L1926
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:37 within `advance_check2`
  store {}* inttoptr (i64 4446596944 to {}*), {}** %.sub, align 8
  %702 = call nonnull {}* @ijl_apply_generic({}* inttoptr (i64 4710351008 to {}*), {}** nonnull %.sub, i32 1)
  call void @ijl_throw({}* %702)
  unreachable

oob:                                              ; preds = %L949
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:26 within `advance_check2`
; ┌ @ array.jl:924 within `getindex`
   %703 = alloca i64, align 8
   store i64 %value_phi24, i64* %703, align 8
   call void @ijl_bounds_error_ints({}* %25, i64* nonnull %703, i64 1)
   unreachable

idxend:                                           ; preds = %L949
   %704 = load double*, double** %374, align 8
   %705 = getelementptr inbounds double, double* %704, i64 %.pre850
   %706 = load double, double* %705, align 8
; └
; ┌ @ operators.jl:591 within `*` @ float.jl:385
   %707 = fmul double %706, %377
; │ @ operators.jl:591 within `*`
   %708 = call nonnull {}* @"j_*_2822"(double %707, {}* nonnull %376) #0
   store {}* %708, {}** %28, align 16
; └
  %709 = call nonnull {}* @"j_+_2823"({}* nonnull %30, {}* nonnull %708) #0
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ abstractarray.jl:95 within `axes`
; ││┌ @ array.jl:151 within `size`
     %710 = load i64, i64* %213, align 8
; │└└
; │┌ @ broadcast.jl:284 within `instantiate`
; ││┌ @ broadcast.jl:543 within `check_broadcast_axes`
; │││┌ @ abstractarray.jl:95 within `axes`
; ││││┌ @ array.jl:151 within `size`
       %711 = bitcast {}* %709 to { i8*, i64, i16, i16, i32 }*
       %712 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %711, i64 0, i32 1
       %713 = load i64, i64* %712, align 8
; │││└└
; │││┌ @ broadcast.jl:540 within `check_broadcast_shape`
; ││││┌ @ broadcast.jl:518 within `_bcsm`
; │││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
        %714 = icmp ne i64 %710, %713
; │││││└
       %715 = icmp ne i64 %713, 1
; ││││└
      %716 = and i1 %714, %715
      br i1 %716, label %L973, label %L992

oob29:                                            ; preds = %L1106
; └└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:27 within `advance_check2`
; ┌ @ array.jl:924 within `getindex`
   %717 = alloca i64, align 8
   store i64 %value_phi24, i64* %717, align 8
   call void @ijl_bounds_error_ints({}* %25, i64* nonnull %717, i64 1)
   unreachable

idxend30:                                         ; preds = %L1106
   %718 = load double*, double** %374, align 8
   %719 = getelementptr inbounds double, double* %718, i64 %.pre850
   %720 = load double, double* %719, align 8
; └
; ┌ @ operators.jl:591 within `*` @ float.jl:385
   %721 = fmul double %720, %377
; │ @ operators.jl:591 within `*`
   %722 = call nonnull {}* @"j_*_2825"(double %721, {}* nonnull %380) #0
   store {}* %722, {}** %28, align 16
; └
  %723 = call nonnull {}* @"j_+_2826"({}* nonnull %36, {}* nonnull %722) #0
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ abstractarray.jl:95 within `axes`
; ││┌ @ array.jl:151 within `size`
     %724 = load i64, i64* %267, align 8
; │└└
; │┌ @ broadcast.jl:284 within `instantiate`
; ││┌ @ broadcast.jl:543 within `check_broadcast_axes`
; │││┌ @ abstractarray.jl:95 within `axes`
; ││││┌ @ array.jl:151 within `size`
       %725 = bitcast {}* %723 to { i8*, i64, i16, i16, i32 }*
       %726 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %725, i64 0, i32 1
       %727 = load i64, i64* %726, align 8
; │││└└
; │││┌ @ broadcast.jl:540 within `check_broadcast_shape`
; ││││┌ @ broadcast.jl:518 within `_bcsm`
; │││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
        %728 = icmp ne i64 %724, %727
; │││││└
       %729 = icmp ne i64 %727, 1
; ││││└
      %730 = and i1 %728, %729
      br i1 %730, label %L1130, label %L1149

oob34:                                            ; preds = %L1263
; └└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:28 within `advance_check2`
; ┌ @ array.jl:924 within `getindex`
   %731 = alloca i64, align 8
   store i64 %value_phi24, i64* %731, align 8
   call void @ijl_bounds_error_ints({}* %25, i64* nonnull %731, i64 1)
   unreachable

idxend35:                                         ; preds = %L1263
   %732 = load double*, double** %374, align 8
   %733 = getelementptr inbounds double, double* %732, i64 %.pre850
   %734 = load double, double* %733, align 8
; └
; ┌ @ operators.jl:591 within `*` @ float.jl:385
   %735 = fmul double %734, %377
; │ @ operators.jl:591 within `*`
   %736 = call nonnull {}* @"j_*_2828"(double %735, {}* nonnull %383) #0
   store {}* %736, {}** %28, align 16
; └
  %737 = call nonnull {}* @"j_+_2829"({}* nonnull %41, {}* nonnull %736) #0
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ abstractarray.jl:95 within `axes`
; ││┌ @ array.jl:151 within `size`
     %738 = load i64, i64* %321, align 8
; │└└
; │┌ @ broadcast.jl:284 within `instantiate`
; ││┌ @ broadcast.jl:543 within `check_broadcast_axes`
; │││┌ @ abstractarray.jl:95 within `axes`
; ││││┌ @ array.jl:151 within `size`
       %739 = bitcast {}* %737 to { i8*, i64, i16, i16, i32 }*
       %740 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %739, i64 0, i32 1
       %741 = load i64, i64* %740, align 8
; │││└└
; │││┌ @ broadcast.jl:540 within `check_broadcast_shape`
; ││││┌ @ broadcast.jl:518 within `_bcsm`
; │││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
        %742 = icmp ne i64 %738, %741
; │││││└
       %743 = icmp ne i64 %741, 1
; ││││└
      %744 = and i1 %742, %743
      br i1 %744, label %L1287, label %L1306

oob39:                                            ; preds = %L1421
; └└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:31 within `advance_check2`
; ┌ @ array.jl:924 within `getindex`
   %745 = alloca i64, align 8
   store i64 %value_phi24, i64* %745, align 8
   call void @ijl_bounds_error_ints({}* %16, i64* nonnull %745, i64 1)
   unreachable

idxend40:                                         ; preds = %L1421
   %746 = load double*, double** %17, align 8
   %747 = getelementptr inbounds double, double* %746, i64 %.pre850
   %748 = load double, double* %747, align 8
; └
; ┌ @ operators.jl:591 within `*` @ float.jl:385
   %749 = fmul double %748, %377
; │ @ operators.jl:591 within `*`
   %750 = call nonnull {}* @"j_*_2831"(double %749, {}* nonnull %376) #0
   store {}* %750, {}** %28, align 16
; └
  %751 = call nonnull {}* @"j_+_2832"({}* nonnull %43, {}* nonnull %750) #0
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ abstractarray.jl:95 within `axes`
; ││┌ @ array.jl:151 within `size`
     %752 = load i64, i64* %45, align 8
; │└└
; │┌ @ broadcast.jl:284 within `instantiate`
; ││┌ @ broadcast.jl:543 within `check_broadcast_axes`
; │││┌ @ abstractarray.jl:95 within `axes`
; ││││┌ @ array.jl:151 within `size`
       %753 = bitcast {}* %751 to { i8*, i64, i16, i16, i32 }*
       %754 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %753, i64 0, i32 1
       %755 = load i64, i64* %754, align 8
; │││└└
; │││┌ @ broadcast.jl:540 within `check_broadcast_shape`
; ││││┌ @ broadcast.jl:518 within `_bcsm`
; │││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
        %756 = icmp ne i64 %752, %755
; │││││└
       %757 = icmp ne i64 %755, 1
; ││││└
      %758 = and i1 %756, %757
      br i1 %758, label %L1446, label %L1465

oob44:                                            ; preds = %L1579
; └└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:32 within `advance_check2`
; ┌ @ array.jl:924 within `getindex`
   %759 = alloca i64, align 8
   store i64 %value_phi24, i64* %759, align 8
   call void @ijl_bounds_error_ints({}* %16, i64* nonnull %759, i64 1)
   unreachable

idxend45:                                         ; preds = %L1579
   %760 = load double*, double** %17, align 8
   %761 = getelementptr inbounds double, double* %760, i64 %.pre850
   %762 = load double, double* %761, align 8
; └
; ┌ @ operators.jl:591 within `*` @ float.jl:385
   %763 = fmul double %762, %377
; │ @ operators.jl:591 within `*`
   %764 = call nonnull {}* @"j_*_2834"(double %763, {}* nonnull %380) #0
   store {}* %764, {}** %28, align 16
; └
  %765 = call nonnull {}* @"j_+_2835"({}* nonnull %99, {}* nonnull %764) #0
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ abstractarray.jl:95 within `axes`
; ││┌ @ array.jl:151 within `size`
     %766 = load i64, i64* %101, align 8
; │└└
; │┌ @ broadcast.jl:284 within `instantiate`
; ││┌ @ broadcast.jl:543 within `check_broadcast_axes`
; │││┌ @ abstractarray.jl:95 within `axes`
; ││││┌ @ array.jl:151 within `size`
       %767 = bitcast {}* %765 to { i8*, i64, i16, i16, i32 }*
       %768 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %767, i64 0, i32 1
       %769 = load i64, i64* %768, align 8
; │││└└
; │││┌ @ broadcast.jl:540 within `check_broadcast_shape`
; ││││┌ @ broadcast.jl:518 within `_bcsm`
; │││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
        %770 = icmp ne i64 %766, %769
; │││││└
       %771 = icmp ne i64 %769, 1
; ││││└
      %772 = and i1 %770, %771
      br i1 %772, label %L1604, label %L1623

oob49:                                            ; preds = %L1737
; └└└└
;  @ /Users/swilliamson/Documents/GitHub/eddy-stresses/checking2.jl:33 within `advance_check2`
; ┌ @ array.jl:924 within `getindex`
   %773 = alloca i64, align 8
   store i64 %value_phi24, i64* %773, align 8
   call void @ijl_bounds_error_ints({}* %16, i64* nonnull %773, i64 1)
   unreachable

idxend50:                                         ; preds = %L1737
   %774 = load double*, double** %17, align 8
   %775 = getelementptr inbounds double, double* %774, i64 %.pre850
   %776 = load double, double* %775, align 8
; └
; ┌ @ operators.jl:591 within `*` @ float.jl:385
   %777 = fmul double %776, %377
; │ @ operators.jl:591 within `*`
   %778 = call nonnull {}* @"j_*_2837"(double %777, {}* nonnull %383) #0
   store {}* %778, {}** %28, align 16
; └
  %779 = call nonnull {}* @"j_+_2838"({}* nonnull %155, {}* nonnull %778) #0
; ┌ @ broadcast.jl:868 within `materialize!` @ broadcast.jl:871
; │┌ @ abstractarray.jl:95 within `axes`
; ││┌ @ array.jl:151 within `size`
     %780 = load i64, i64* %157, align 8
; │└└
; │┌ @ broadcast.jl:284 within `instantiate`
; ││┌ @ broadcast.jl:543 within `check_broadcast_axes`
; │││┌ @ abstractarray.jl:95 within `axes`
; ││││┌ @ array.jl:151 within `size`
       %781 = bitcast {}* %779 to { i8*, i64, i16, i16, i32 }*
       %782 = getelementptr inbounds { i8*, i64, i16, i16, i32 }, { i8*, i64, i16, i16, i32 }* %781, i64 0, i32 1
       %783 = load i64, i64* %782, align 8
; │││└└
; │││┌ @ broadcast.jl:540 within `check_broadcast_shape`
; ││││┌ @ broadcast.jl:518 within `_bcsm`
; │││││┌ @ range.jl:1111 within `==` @ promotion.jl:477
        %784 = icmp ne i64 %780, %783
; │││││└
       %785 = icmp ne i64 %783, 1
; ││││└
      %786 = and i1 %784, %785
      br i1 %786, label %L1762, label %L1781
; └└└└
}
