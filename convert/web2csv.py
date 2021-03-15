import pandas as pd


web_f = ''
web_f_csv = ''
columns = ['dens', 'Vx', 'Vy', 'Vz', 'l1', 'l2', 'l3', 'e1x', 'e1y', 'e1z', 'e2x', 'e2y', 'e2z', 'e3x', 'e3y', 'e3z']

vweb = pd.read_csv(web_f, separator='\t')

vweb.to_csv()

