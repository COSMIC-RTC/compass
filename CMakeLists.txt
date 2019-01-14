cmake_minimum_required (VERSION 3.8)

project(compass
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS library"
    LANGUAGES CUDA CXX
)

# WYRM IMPORT #
################

if(NOT TARGET wyrm)
	if(NOT DEFINED WYRM_ROOT)
		set(WYRM_ROOT $ENV{WYRM_ROOT})
	endif()
	if(IS_DIRECTORY ${WYRM_ROOT})
		add_subdirectory(${WYRM_ROOT} wyrm_dir)
    else()
        message(ERROR "wyrm path not found, specify it by set WYRM_ROOT to the right path" )
    endif()
endif()

add_subdirectory(libcarma)
add_subdirectory(libsutra)
# add_subdirectory(shesha)
