import pymysql.cursors

class db_conn(object):

	def __init__(self,database='',server=None,test=False):

		self.readOnly = 0

		#define database connection info
		self.db = database
		if database == 'koa':
			self.server			= 'koaserver.keck.hawaii.edu'
			self.backupServer	= ''
			self.user			= 'koaup'
			self.pwd			= 'kexpress_ru'

			if test==True:
				self.db = 'test'
		else:
			return

	def db_connect(self):
		'''
		Connect to the specified database.  If primary server is down and backup
		is specified, then connect to it.  This also set the readOnly flag to 1.
		'''

		self.dbConn = pymysql.connect(self.server, self.user, self.pwd, self.db, cursorclass=pymysql.cursors.DictCursor)

	def db_close(self):
		'''
		Closes the current database connection
		'''

		if self.dbConn:
			self.dbConn.close()

	def do_query(self, query, output):
		self.db_connect()

		# Save as a list of dictionaries
		query = ''.join(query)
		result = ''
		with self.dbConn.cursor() as cursor:
			num = cursor.execute(query)
			result = cursor.fetchall()
			# if num == 0:
			# 	value = {}
			# 	if output == 'txt':
			# 		value = ' '
			# else:
			# 	value = []
			# 	count = 0
			# 	for row in result:
			# 		value.insert(count, {})
			# 		for f in row.keys():
			# 			value[count][f] = str(row[f])
			# 		count += 1

			# 	# Convert to text output
			# 	if output == 'txt':
			# 		v = ''
			# 		count = 0
			# 		for row in value:
			# 			if count > 0:
			# 				v = v + '<br>'
			# 			for key, val in row.items():
			# 				v = v + ' ' + str(val)
			# 			count += 1
			# 		value = v

		self.db_close()
		return result
