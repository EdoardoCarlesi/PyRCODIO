
'''
	(file_lg_name, file_sub_name) = settings.get_zoom_output()

	file_lg_txt = open(file_lg_name, 'wb')
	file_lg_txt.write(file_lg_header)
	file_sub_txt = open(file_sub_name, 'wb')
	file_sub_txt.write(file_sub_header)

	# Subhalo mass function variables
	file_png_mfs1 = settings.get_png_output('mf_01')
	file_png_mfs2 = settings.get_png_output('mf_02')
	x_mf1 = []; 	y_mf1 = []
	x_mf2 = []; 	y_mf2 = []
	n_mf = 0
		file_png_name = settings.get_png_output('lg')
		this_file_ahf = settings.ahf_path + ahf_snap; 		this_file_ahf_alt = settings.ahf_path + ahf_snap_alt
		this_file_gad = settings.file_z0
		print_run = base_run + '_' + run_num
		good_lgs = 0
		
		if os.path.exists(this_file_ahf):
				print 'Reading in AHF file: ', this_file_ahf
				ahf_all = read_ahf(this_file_ahf)
				these_lg = find_lg(ahf_all, lg_model)
				n_lgs = int(len(these_lg))
		else:
				print 'No AHF file: ', this_file_ahf
				n_lgs = 0;
		
		print 'Found a total of %d LG pairs.' % (n_lgs)

		if n_lgs > 0:
			# If there are more candidates we need to find the right one
			rating = 1000
			for ind in range(0, n_lgs):
				lg = these_lg[ind]
				lg.c_box = settings.box_center

				if lg.rating() < rating and (lg.LG1.npart > part_min) and (lg.LG2.npart > part_min):
					print 'Old rating: %f new rating %f this index %d' % (rating, lg.rating, ind)
					good_lgs += 1
					rating = lg.rating
					best_lg = lg

			if good_lgs > 0:
				print 'Best LG: ', best_lg.info()

				# Pair properties
				file_lg_line = best_lg.info()

				#if vel_r < vrad_max and m12 < m_max:
				file_lg_txt.write(file_lg_line)

				# Now take care of the substructure
				these_sub1 = find_halos(best_lg.LG1, ahf_all, fac_r * best_lg.LG1.r)
				subs1 = SubHalos(best_lg.LG1, these_sub1, print_run, 'M31')
				subs1.anisotropy("part", np_sub_min)
				subs1.basis_eigenvectors("inertia")

				these_sub2 = find_halos(best_lg.LG2, ahf_all, fac_r * best_lg.LG2.r)
				subs2 = SubHalos(best_lg.LG2, these_sub2, print_run, 'MW')
				subs2.anisotropy("part", np_sub_min)
				subs2.basis_eigenvectors("inertia")

				file_sub_line  = subs1.all_info("part", 2 * np_sub_min)
				file_sub_line += subs2.all_info("part", 2 * np_sub_min)
				file_sub_line += '\n ------------------------ \n'
				file_sub_txt.write(file_sub_line)

				(x_m, y_n) = subs1.mass_function()
				x_mf1.append(x_m); 	y_mf1.append(y_n)

				(x_m, y_n) = subs2.mass_function()
				x_mf2.append(x_m); 	y_mf2.append(y_n)
				n_mf += 1
			
			# Do some plots at the end
			if do_plots == "true":
				plot_lg(this_file_gad, file_png_name, best_lg.LG1, best_lg.LG2, reduce_fac, 1, plot_pos)
	
				if (subrun_i == subrun_end-1):
					plot_massfunctions(x_mf1, y_mf1, n_mf, file_png_mfs1)
					plot_massfunctions(x_mf2, y_mf2, n_mf, file_png_mfs2)
					del x_m	;	del y_n
					del x_mf1; 	del x_mf2
					del y_mf1; 	del y_mf2
					n_mf = 0

			# Compute halo evolution
'''
			

