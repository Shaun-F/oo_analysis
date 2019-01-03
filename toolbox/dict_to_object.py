class dict_to_object:
	def __init__(self, **vars):
		self.__dict__.update(vars)
		