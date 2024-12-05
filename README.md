#1 Run RosettaMatch to search for potential metal coordination sites 
./path/RosettaMatch -s input_scaffold/test_scaffold.pdb -match:scaffold_active_site_residues_for_geomcsts params_file/test_scaffold.pos -match::lig_name ZNX -match::geometric_constraint_file params_file/ZNX_HD.cst -match::output_matchres_only -extra_res_fa params_file/ZNX.params -output_matches_per_group 1 -packing:ex1 -packing:ex2 -euclid_bin_size 1 -euler_bin_size 10.0 -bump_tolerance 0.5 -match:consolidate_matches -match:output_matchres_only -match:output_format PDB;

#2 List RosettaMatch outputs
python script/List_match_outputs.py -input_list test_scaffold/list -output_path test_scaffold -output_name test_scaffold

#3 List coordination frames
python script/List_coordination_frames.py -scaffold_pdb input_scaffold/test_scaffold.pdb -matchlist test_scaffold/test_scaffold.list -output_path test_scaffold -output_name test_scaffold

#4 Protein Dock
./path/interface_metal_homooligomer -use_ss -motif_pos 1 -num_repeats 1 -rpx_db params_file/rpx_cart2.0_ang26.0_ss_ALLAA.phmap -symmetry C3 -output_dir test_scaffold_dock_output/ -len 88 -ss HHHHHHHHHHHHHHHHHHHLLHHHHHHHHHHHHHHHHHHHHLLHHHHHHHHHHHHHHHHHHHLLLHHHHHHHHHHHHHHHHHHHHHHH -interface_metal_scaffold_pdb input_scaffold/test_scaffold.pdb -output_prefix test__C3 -nstruct 5 -total_score_cutoff 20.0 -inter_chain_score_cutoff -0.5 -sidechain_neighbor_cutoff 0.1 -sampling_stage_cutoff_tolerance_factor 0.1 -interface_metal_score_cutoff 1 -interface_metal_distance_optimization_weight 1 -interface_metal_type zinc -num_interface_metals 1   -interface_metal_config1 test_scaffold.dat  -mcmc_inner_cycles 10000 -mcmc_outer_cycles 3;

#5 Replace metal coordination
python script/Replace_zinc.py -pdb test_scaffold_dock_output/test__C3_1732945764539097525_SN_0.631_TotalSc_-1.16_RPX_-2.13_Clash_0_MetalDist_0.973_Zn_Cluster0127_Zn_Cluster0053.pdb -scaffold_pdb input_scaffold/test_scaffold.pdb -metal_pdb params_file/zn.pdb -metal_params params_file/ZNX.params -match_path test_scaffold/ -coordination_frames test_scaffold/test_scaffold.dat -output_path replace_zinc/ 

#6 Generate cst file
python script/constrain_generation.py -pdb replace_zinc/test__C3_1732945764539097525_SN_0.631_TotalSc_-1.16_RPX_-2.13_Clash_0_MetalDist_0.973_Zn_Cluster0127_Zn_Cluster0053.pdb -coordination_frame test_scaffold/test_scaffold.dat -scaffold_pdb input_scaffold/test_scaffold.pdb -metal_num 1  -output_path constrain_and_asym_pdb -symmetry_number 3 -match_path test_scaffold

#7 Find designable residues
python script/find_designable_res.py -pdb replace_zinc/test__C3_1732945764539097525_SN_0.631_TotalSc_-1.16_RPX_-2.13_Clash_0_MetalDist_0.973_Zn_Cluster0127_Zn_Cluster0053.pdb -not_design_threshold 8 -full_pose -output constrain_and_asym_pdb/designable.list -symmetry C3 -num_zincs 1 -pos_file params_file/test_scaffold.pos

#8 Fastdesign
./path/Rosetta_scripts -s constrain_and_asym_pdb/test__C3_1732945764539097525_SN_0.631_TotalSc_-1.16_RPX_-2.13_Clash_0_MetalDist_0.973_Zn_Cluster0127_Zn_Cluster0053.pdb  -parser:protocol Fastdesign.xml -out:file:renumber_pdb true -corrections:beta_nov16 -parser:script_vars sym_file=params_file/C3_Z.sym  cst_file=constrain_and_asym_pdb/test__C3_1732945764539097525_SN_0.631_TotalSc_-1.16_RPX_-2.13_Clash_0_MetalDist_0.973_Zn_Cluster0127_Zn_Cluster0053.cst designable_res=23,45,49,53,67,68,70,71,74,75,78,79,81,82,83,85,112,134,138,142,156,157,159,160,163,164,167,168,170,171,172,174,201,223,227,231,245,246,248,249,252,253,256,257,259,260,261,263 -output_pose_energies_table false -mute protocols.rosetta_scripts.ParsedProtocol.REPORT -extra_res_fa params_file/ZNX_virtual.params -out:pdb_gz  -output_virtual -nstruct 1;
