

def gen_table_lg():
	sql_table = """
	CREATE TABLE local_group (
	lg_number INTEGER PRIMARY KEY,
	simu_code VARCHAR(10),
	IDMW INT64;
	IDM31 INT64;
	v_rad DOUBLE(10,4),
	r_mwm31 DOUBLE(10,4),
	m_ratio DOUBLE(10,4),
	VirgoM DOUBLE(10,4),
	VirgoX DOUBLE(10,4),
	VirgoY DOUBLE(10,4),
	VirgoZ DOUBLE(10,4),
	M31_M DOUBLE(10,4),
	M31_X DOUBLE(10,4),
	M31_Y DOUBLE(10,4),
	M31_Z DOUBLE(10,4),
	MW_M DOUBLE(10,4),
	MW_X DOUBLE(10,4),
	MW_Y DOUBLE(10,4),
	MW_Z DOUBLE(10,4),
	COM_X DOUBLE(10,4),
	COM_Y DOUBLE(10,4),
	COM_Z DOUBLE(10,4)
	);
	"""


