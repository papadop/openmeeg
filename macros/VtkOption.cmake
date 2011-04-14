#------------------------------------------------------------
# VTK library
#------------------------------------------------------------

OPTION(USE_VTK "Build the project using VTK" OFF)
MARK_AS_ADVANCED(USE_VTK)

IF (USE_VTK)
    FIND_PACKAGE(VTK)
    IF (VTK_FOUND)
        SET(USE_VTK 1)
        INCLUDE(${VTK_USE_FILE})
        SET (VTK_LIBRARIES
            # vtkRendering
            vtkGraphics
            # vtkHybrid
            # vtkImaging
            vtkIO
            # vtkFiltering
            # vtkGenericFiltering
            vtkCommon
            # vtkDICOMParser
            # vtkzlib
        )
        SET(OPENMEEG_LIBRARY_DIRS 
            ${OPENMEEG_LIBRARY_DIRS}
            ${VTK_LIBRARY_DIRS}
        )
        SET(OPENMEEG_INCLUDE_DIRECTORIES
            ${OPENMEEG_INCLUDE_DIRECTORIES}
            ${VTK_INCLUDE_DIRS}
        )
        SET(OPENMEEG_LIBRARIES
            ${OPENMEEG_LIBRARIES}
            ${VTK_LIBRARIES}
        )
    ELSE()
        MESSAGE(FATAL_ERROR "Please set VTK_DIR")
    ENDIF()
ENDIF()
