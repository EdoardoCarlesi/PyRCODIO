import pandas as pd
import sqlite3
import sys
from .mtree import *

'''
	SQL_IO works as an interface between a SQL database containing tables of halos and a python program
'''

class SQL_IO:
	dbName = ''
	nSnaps = 0
	
	# Initialize as dummy objects, just to have the right type
	mydb = sqlite3.connect('')
	cursor = sqlite3.connect('').cursor()

	# Pandas data frame 
	dataframe = pd.DataFrame()

	def __init__(self, dbName, nSnaps):	
		self.dbName = dbName	
		self.nSnaps = nSnaps
		self.mydb = sqlite3.connect(self.dbName)
		self.cursor = self.mydb.cursor()

		print('Reading DB into pandas dataframe...')
		selectStr = "SELECT * FROM halo"
		self.dataframe = pd.read_sql_query(selectStr, self.mydb)
		print('Done.')
	
	# Create a halo table into the database
	def halo_table(self):
		'Create halo table for SQL database'

		createStr = """CREATE TABLE IF NOT EXISTS halo (
				haloID INT64,
				simuCode VARCHAR(10), 
				allNumPart INT ARRAY[""" + str(self.nSnaps) + """], 
				allHaloIDs INT64 ARRAY[""" + str(self.nSnaps) + """],
				X FLOAT ARRAY[""" + str(self.nSnaps) + """], 
				Y FLOAT ARRAY[""" + str(self.nSnaps) + """], 
				Z FLOAT ARRAY[""" + str(self.nSnaps) + """], 
				VX FLOAT ARRAY[""" + str(self.nSnaps) + """], 
				VY FLOAT ARRAY[""" + str(self.nSnaps) + """], 
				VZ FLOAT ARRAY[""" + str(self.nSnaps) + """]
			)"""

		self.cursor.execute(createStr)
	
	# Update the properties of a given halo inside a database
	def update_tree(self, ID):
		strCheck = self.select_tree(ID)
		
		if strCheck != None:
			print('Updating the database %s for HaloID %s.' % (self.dbName, ID) )
		
		'TODO: add some halo (or progenitor) properties: positions, velocties, etc.'
	
	# Remove an halo from the database
	def delete_tree(self, ID):
		deleteStr = """DELETE FROM halo WHERE haloID=""" + str(ID)
		self.cursor.execute(deleteStr)

	'''
	DONT write stuff into the database from PyRCODIO
	# Add a new tree to the database
	def insert_tree(self, ID, simuCode, haloParts, haloIDs):
		checkStr = None

		haloPartsStr = ''
		haloIDsStr = ''
	
		for iPart in range (0, len(haloParts)):
			if iPart < len(haloParts) - 1:
				haloPartsStr += str(haloParts[iPart]) + ', '
				haloIDsStr += str(haloIDs[iPart]) + ', '
			else:	
				haloPartsStr += str(haloParts[iPart])
				haloIDsStr += str(haloIDs[iPart])

		insertStr = """INSERT INTO halo (haloID, allNumPart, allHaloIDs, simuCode) 
				VALUES ('"""+ ID +"""', '""" + simuCode + """',  
				'""" + haloPartsStr + """', '""" + haloIDsStr + """' );"""
		self.cursor.execute(insertStr)
	'''

	# Simple database query
	def select_tree(self, ID, columnName):
		try:
			thisStr = self.dataframe.loc[self.dataframe['haloID'] == ID, columnName].values[0].split(", ")
			nStr = len(thisStr)

			if columnName == 'allHaloIDs':
				thisTree = np.zeros((nStr), dtype=np.ulonglong)
			elif columnName == 'allNumPart':
				thisTree = np.zeros((nStr), dtype=np.uint)

			iTree = 0

			for iStr in thisStr:
				iStr.replace('u', '')
				thisTree[iTree] = long(iStr)
				#print iTree, iStr
				iTree += 1

			return thisTree
		except:
			#print(ID, ' not found')
			return np.zeros((1))


	# Simple database query
	def select_trees(self, IDs):
		allTrees = []

		iID = 0
		for ID in IDs:
			thisTree = self.dataframe.loc[self.dataframe['haloID'] == ID, 'allHaloIDs'].values[0].split(", ")
			allTrees.append(thisTree)

		return allTrees


	# Given a halo ID, retrieve the full halo merger history from the database
	def get_full_mtree(self, ID):
		mergerTree = MTree(self.nSnaps, ID)
		
		tree_line = self.select_tree(ID)
		these_ids = tree_line[1].split()		
		these_nps = tree_line[3].split()		
		ids = []; nps = []; nPt = 0

		# We loop on 54 values (which is what we expect) though for smaller haloes the tree might be shorter
		for nPt in range(0, self.nSnaps):
			
			try:
				this_id = these_ids[nPt].replace(",", "")
				this_np = these_nps[nPt].replace(",", "")
			except:
				'This tree breaks up before expected. Maybe it is a small halo, so lets just add token values to it.'
	
			if this_np == '' or this_np == ' ':
				this_np = 0	

			if this_id == '' or this_id == ' ':
				this_id = ids[nPt-1]
	
			ids.append(this_id)
			nps.append(this_np)

		mergerTree.fill_mass_id(nps, ids)

		return mergerTree

	# Close the database
	def close(self):
		self.mydb.close()



