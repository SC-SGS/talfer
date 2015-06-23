/*
 * CVTKTools.hpp
 *
 *  Created on: 10.04.2015
 *      Author: axel
 */

#ifndef CVTKTOOLS_HPP_
#define CVTKTOOLS_HPP_

#include <iostream>
#include <stdio.h>
#include <fstream>
#include "CDataArray2D.hpp"


static void writeVTKFile(CDataArray2D<float, 9> visData, int timestep)
{
	char numstr[21];
	sprintf(numstr, "%d", timestep);

	int width = visData.width;
	int height = visData.height;

	printf("%i, %i", width, height);

	std::string filename;
	filename.append("./");
	filename.append("output");
	filename.append("/");
	filename.append("LBMGrid_");
	filename.append(numstr);
	filename.append(".vtr");
		// genearte filebuf to write to file
	std::filebuf fb;
	fb.open(const_cast<char *>(filename.c_str()), std::ios::out);
	std::ostream os(&fb);

	if(!fb.is_open()){
		printf("error");
	}else{
		//printf("%s geöffnet", const_cast<char *>(filename.c_str()));
	}

	// vtk Header
	os << "<?xml version=\"1.0\"?>" << std::endl;
	os << "<VTKFile type=\"RectilinearGrid\">" << std::endl;

	// domain size
	os << "<RectilinearGrid WholeExtent=\"0 " << width << " 0 " << height
	   << " 0 0\">" << std::endl;
	os << "<Piece Extent=\"" << 0 << " " << width << " " << 0 << " " << height << " 0 0\">" << std::endl;

	// specify coordinates
	os << "<Coordinates>" << std::endl;
	os << "<DataArray type=\"Float64\" format=\"ascii\">" << std::endl;

	for(int i = 0; i <= width; ++i){
		os << (float)i << " ";
	}
	os << std::endl;

	os << "</DataArray>" << std::endl;

	// along y-Axis
	os << "<DataArray type=\"Float64\" format=\"ascii\">" << std::endl;
	for(int i = 0; i <= height; ++i){
		os << (float)i << " ";
	}
	os << std::endl;
	os << "</DataArray>" << std::endl;

	// along z-Axis (nothing to do => 2D)
	os << "<DataArray type=\"Float64\" format=\"ascii\">0 0</DataArray>"
	   << std::endl;

	os << "</Coordinates>" << std::endl;

	// add payload
	os << "<CellData Vectors=\"velocity\" Scalars=\"rho, p_force\">"
	   << std::endl;

	// 1 velocity field
	os << "<DataArray Name=\"velocity\" NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">"
	  << std::endl;

	for(int j = height-1; j >=0; --j){
	    for(int i = 0; i < width; ++i){
	    	//printf("tut was bei %i, %i", i, j);
			os << std::scientific << visData.getRef(i, j, 0) << " " << visData.getRef(i, j, 1) << " " << 0.  // war mal u[i][j] v[i][j]
			   << std::endl;
		}
	}

	os << "</DataArray>" << std::endl;

	// 2 pressure
	os << "<DataArray type=\"Float64\" Name=\"rho\" format=\"ascii\">"
	   << std::endl;

	for(int j = height-1; j >=0; --j){
		for(int i = 0; i < width; ++i){
			os << std::scientific << visData.getRef(i, j, 2) << std::endl;
		}
	}

	os << "</DataArray>" << std::endl;

	// 3 vorticity
	os << "<DataArray type=\"Float64\" Name=\"p_force\" format=\"ascii\">"
	   << std::endl;

	for(int j = height-1; j >=0; --j){
		for(int i = 0; i < width; ++i){
			os << std::scientific << visData.getRef(i, j, 3) << std::endl;
		}
	}

	os << "</DataArray>" << std::endl;

	os << "</CellData>" << std::endl;

	os << "</Piece>" << std::endl;

	os << "</RectilinearGrid>" << std::endl;
	os << "</VTKFile>" << std::endl;

}


#endif /* CVTKTOOLS_HPP_ */
