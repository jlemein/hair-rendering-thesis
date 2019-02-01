/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   OpenVdbReader.cpp
 * Author: jeffrey
 *
 * Created on February 1, 2019, 4:25 PM
 */

#include <openvdb/openvdb.h>

#include "OpenVdbReader.h"
#include "../hairstruct.h"

OpenVdbReader::OpenVdbReader(std::string fileName) : mInputFileName(fileName) {
}

OpenVdbReader::OpenVdbReader(const OpenVdbReader& orig) {
}

OpenVdbReader::~OpenVdbReader() {
}

void OpenVdbReader::initialize() {
    openvdb::initialize();
    openvdb::io::File file(mInputFileName);
    file.open();

    for (openvdb::io::File::NameIterator nameIter = file.beginName();
            nameIter != file.endName(); ++nameIter) {
        std::cout << "Found grid with name: " << nameIter.gridName() << std::endl;
        openvdb::GridBase::Ptr grid = file.readGrid(nameIter.gridName());
        printMetaDataForGrid(grid);

        mGrids.push_back(grid);

        if (nameIter.gridName() == "HairDensityVolume") {
            mHairDensityGrid = grid;
        }
    }

    file.close();
}

/**
 * Interpolates the voxel grid from a 3D point in space to another point
 * in 3D space
 * @param from
 * @param to
 * @return
 */
float OpenVdbReader::interpolate(const Point3& from, const Point3 to) {
    std::cout << "[WARNING]: interpolate(Point3&, Point$) is not implemented\n";
    return -1.0f;
}

float OpenVdbReader::interpolateToInfinity(const Point3& from, const Point3 direction) {
    std::cout << "[WARNING]: interpolateToInfinity(Point3&, Point$) is not implemented\n";
    return -1.0f;
}

void OpenVdbReader::printMetaDataForGrid(openvdb::GridBase::Ptr grid) const {
    std::cout << "Metadata for grid with name '" << grid->getName() << "':\n";
    for (openvdb::MetaMap::MetaIterator iter = grid->beginMeta();
            iter != grid->endMeta(); ++iter) {
        const std::string& name = iter->first;
        openvdb::Metadata::Ptr value = iter->second;
        std::string valueAsString = value->str();
        std::cout << "\t" << name << " = " << valueAsString << std::endl;
    }
}

void OpenVdbReader::printMetaDataForAllGrids() const {
    for (auto& grid : mGrids) {
        printMetaDataForGrid(grid);
        std::cout << std::endl;
    }
}

void OpenVdbReader::printMetaDataForHairDensityGrid() const {
    if (mHairDensityGrid != 0) {
        printMetaDataForGrid(mHairDensityGrid);
    } else {
        std::cout << "Hair density grid does not exist\n";
    }
}

void OpenVdbReader::printVdbInfo() const {
    for (auto& grid : mGrids) {
        std::cout << "Voxel info for grid with name '" << grid->getName() << "'\n\n";
        printMetaDataForGrid(grid);
    }
}