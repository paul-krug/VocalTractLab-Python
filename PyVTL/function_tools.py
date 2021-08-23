import os
import itertools
import warnings

#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
# 		input argument related functions
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def check_if_input_lists_are_valid( input_lists, instances_list ):
	valid_lists = []
	for input_list, instance in zip( input_lists, instances_list ):
		valid_lists.append( check_if_list_is_valid( input_list, instance ) )
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
def check_if_list_is_valid( input_list, instance ):
	is_iterable = check_if_object_is_iterable( input_list )
	if isinstance( input_list, str ) or is_iterable == False:
		warnings.warn( 'input is either not iterable or a single string. The input gets turned into a list now.' )
		input_list = [ input_list ]
	if input_list and all( isinstance( x, instance ) for x in input_list ):
		return input_list
	else:
		raise TypeError( 'a list containing a non-{} type object was passed, but list of {} was expected.'.format( instance, instance ) )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def check_if_all_elements_are_equal( iterable ):
	g = itertools.groupby(iterable)
	return next(g, True) and not next(g, False)
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def make_output_path( query_path, output_path ):
	if '.' not in file_path:
		raise ValueError( 'The path: {} is not valid because it does not specify the file extension!'.format( output_path ) )
	if query_path in (None, ''):
		file_extension = file_path.rsplit('.')[1]
		index = 0
		while os.path.exists( output_path ):
			index += 1
			output_path = output_path.replace( file_extension, '_{}.{}'.format( index, file_extension ) )
		log.info( 'No output file path for audio file was specified, saving file to {}'.format( audio_file_path ) )
	if not os.path.exists( os.path.dirname( output_path ) ) and  os.path.dirname( output_path ) not in ( '', ' ', None ):
		os.mkdir( os.path.dirname( output_path ) )
		log.info( 'Output directory {} did not exist and was created.'.format( os.path.dirname( output_path ) ) )
	return output_path
#---------------------------------------------------------------------------------------------------------------------------------------------------#