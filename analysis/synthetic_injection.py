#Scrip injects synthetic axion signals into data
import
import sys 
sys.path.insert(0, os.path.abspath("../")
sys.path.append("../signals")

from signal_lib import signal as sig_lib



class injector():
	
	
	def __init__(self, **kwargs):
		
		#set class definitions
		for name, arg in kwargs:
			setattr(self, name, arg)
			
		self.signal_class = sig_lib()
		
	def generate_signals(self):
		signal_name = self.signal
		signal_shape = getattr(self.signal_class, signal_name)(