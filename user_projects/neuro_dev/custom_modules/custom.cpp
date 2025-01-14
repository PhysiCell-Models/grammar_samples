/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

void create_cell_types(void)
{
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/*
       Cell rule definitions 
	*/

	setup_cell_rules();

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function;

	find_cell_definition("apical")->is_movable = false;

	find_cell_definition( "rgc" )->functions.update_phenotype = rgc_phenotype_function;
	find_cell_definition( "layer_6" )->functions.update_phenotype = migrating_phenotype_function;
	find_cell_definition( "layer_5" )->functions.update_phenotype = migrating_phenotype_function;
	find_cell_definition( "layer_4" )->functions.update_phenotype = migrating_phenotype_function;
	find_cell_definition( "layer_3" )->functions.update_phenotype = migrating_phenotype_function;
	find_cell_definition( "layer_2" )->functions.update_phenotype = migrating_phenotype_function;
	
	find_cell_definition( "layer_6" )->functions.custom_cell_rule = custom_function;
	find_cell_definition( "layer_5" )->functions.custom_cell_rule = custom_function;
	find_cell_definition( "layer_4" )->functions.custom_cell_rule = custom_function;
	find_cell_definition( "layer_3" )->functions.custom_cell_rule = custom_function;
	find_cell_definition( "layer_2" )->functions.custom_cell_rule = custom_function;
	find_cell_definition( "pial" )->functions.custom_cell_rule = pial_function;

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment(void)
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM

	initialize_microenvironment();

	return; 
}

void setup_tissue(void)
{
	setup_tissue_domain();
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();
	for ( int i = 0; i < (*all_cells).size(); i++ )
	{
		Cell* pC = (*all_cells)[i];
		if (pC->type_name == "rgc")
		{
			pC->phenotype.cycle.data.elapsed_time_in_phase = UniformRandom() / get_single_base_behavior(pC, "cycle entry");
		}
	}

	return;
}

void setup_tissue_domain(void)
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{
	// check if any neighbors are pial cells
	// if so, stop migration
	// if not, continue migration
	for (int i = 0; i < pCell->state.neighbors.size(); i++)
	{
		if (pCell->state.neighbors[i]->type_name == "pial")
		{
			phenotype.motility.migration_speed = 0;
			pCell->functions.custom_cell_rule = NULL;
			return;
		}
	}
	return;
}

void rgc_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	return;
}

void migrating_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	if( phenotype.motility.migration_speed > 0.2 )
	{
		phenotype.mechanics.cell_cell_repulsion_strength = 0;
	}
	else
	{
		phenotype.mechanics.cell_cell_repulsion_strength = find_cell_definition( pCell->type_name )->phenotype.mechanics.cell_cell_repulsion_strength; 
	}
	return;

}

// DZ: 6-19 - a pial cell phenotype function to prevent the cells from being caught at the sides
void pial_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	/*
	std::vector<double> boundaries = pCell->get_microenvironment()->mesh.bounding_box;
	float diameter = 2 * pCell->phenotype.geometry.radius;
	if( pCell->position[0] < boundaries[0] + diameter)
	{
		pCell->position[0] = boundaries[0] + diameter;
		// std::cout << "X position updated to " << pCell->position[0] << std::endl;
	}
	else if( pCell->position[0] > boundaries[3] - diameter)
	{
		pCell->position[0] = boundaries[3] - diameter;
		// std::cout << "X position updated to " << pCell->position[0] << std::endl;
	}
	*/
	return;
}

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; }

void asymmetric_division_function(Cell *pC1, Cell *pC2)
{
	double ct = PhysiCell_globals.current_time;
	std::string new_cell_type;
	if (ct <= 1440.0)
	{
		return; // symmetric divsion into two rgc cells
	}
	else if (ct <= parameters.doubles("layer_6_end_time"))
	{
		new_cell_type = "layer_6";
	}
	else if (ct <= parameters.doubles("layer_5_end_time"))
	{
		new_cell_type = "layer_5";
	}
	else if (ct <= parameters.doubles("layer_4_end_time"))
	{
		new_cell_type = "layer_4";
	}
	else if (ct <= parameters.doubles("layer_3_end_time"))
	{
		new_cell_type = "layer_3";
	}
	else if (ct <= 12960)
	{
		new_cell_type = "layer_2";
	}
	else
	{
		return; // no more asymmetric division
	}
	Cell_Definition* pCD = cell_definitions_by_name[new_cell_type];
	pC2->convert_to_cell_definition( *pCD );
	return;
}