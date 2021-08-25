import os
import itertools
import warnings
import logging

logging.basicConfig()
log = logging.getLogger(__name__)
log.setLevel(logging.WARNING)

#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
# 		input argument related functions
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def check_if_input_lists_are_valid( input_lists, instances_list ):
	valid_lists = []
	for input_list, instances in zip( input_lists, instances_list ):
		valid_lists.append( check_if_list_is_valid( input_list, instances ) )
	return valid_lists
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def check_lengths_of_input_lists( input_lists ):
	list_lengths = [ len( input_list ) for input_list in input_lists ]
	if not check_if_all_elements_are_equal( list_lengths ):
		warnings.warn( 'input list do not have the same lengths, shorter lists will be padded with "None".' )
		max_length = max( list_lengths )
		#print( max_length )
		for input_list in input_lists:
			while len( input_list ) < max_length:
				input_list.append( None )
		#print( input_list )
	return input_lists
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def check_if_object_is_iterable( query ):
	try:
		iter( query )
	except TypeError:
		return False
	return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def check_if_list_is_valid( input_list, instances ):
	is_iterable = check_if_object_is_iterable( input_list )
	if isinstance( input_list, str ) or is_iterable == False:
		warnings.warn( 'input is either not iterable or a single string. The input gets turned into a list now.' )
		input_list = [ input_list ]
	if input_list and all( isinstance( x, instances ) for x in input_list ):
		return input_list
	else:
		raise TypeError( 'a list containing a non-{} type object was passed, but list of {} was expected.'.format( instances, instances ) )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def check_if_all_elements_are_equal( iterable ):
	g = itertools.groupby(iterable)
	return next(g, True) and not next(g, False)
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def make_output_path( query_path, output_path ):
	if '.' not in output_path:
		raise ValueError( 'The output path: {} is not valid because it does not specify the file extension!'.format( output_path ) )
	if query_path in (None, ''):
		query_path = output_path
	else:
		output_path = query_path
	file_extension = '.' + query_path.rsplit('.')[1]
	index = 0
	while os.path.exists( output_path ):
		index += 1
		output_path = output_path.replace( file_extension, '_{}{}'.format( index, file_extension ) )
	log.info( 'Saving file to {}'.format( output_path ) )
	if not os.path.exists( os.path.dirname( output_path ) ) and  os.path.dirname( output_path ) not in ( '', ' ', None ):
		os.mkdir( os.path.dirname( output_path ) )
		log.info( 'Output directory {} did not exist and was created.'.format( os.path.dirname( output_path ) ) )
	return output_path
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def make_output_dir( query_dir, output_dir ):
	if '.' in output_dir:
		raise ValueError( 'The output dir: {} is not valid because it contains a dot!'.format( output_dir ) )
	if query_dir in ( None, '', ' ' ):
		query_dir = output_dir
	else:
		output_dir = query_dir
	index = 0
	output_dir_index = output_dir + '_0'
	while os.path.exists( output_dir ):
		index += 1
		output_dir = output_dir_index.replace( '_0', '_{}'.format( index ) )
	log.info( 'Saving file to {}'.format( output_dir ) )
	if not os.path.exists( output_dir ) :
		os.mkdir( output_dir )
		log.info( 'Output directory {} did not exist and was created.'.format( output_dir ) )
	return output_dir
#---------------------------------------------------------------------------------------------------------------------------------------------------#