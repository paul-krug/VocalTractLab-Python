def worker( args ):
	func, arg = args
	func( arg )
	return