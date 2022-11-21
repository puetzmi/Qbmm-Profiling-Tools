get_filename_component(SRC_DIR ${SRC} DIRECTORY)
# Generate a git-describe version string from Git repository tags
execute_process(
COMMAND ${GIT_EXECUTABLE} describe --dirty
WORKING_DIRECTORY ${SRC_DIR}
OUTPUT_VARIABLE GIT_DESCRIBE_VERSION
RESULT_VARIABLE GIT_DESCRIBE_ERROR_CODE
OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT GIT_DESCRIBE_ERROR_CODE)
    set(CMAKE_GIT_BUILD_ID "${COMPILER_INFO}--${GIT_DESCRIBE_VERSION}")
else()
    set(CMAKE_GIT_BUILD_ID "${COMPILER_INFO}--v0.0.0-unknown")
    message(WARNING "Failed to determine CMAKE_GIT_BUILD_ID from Git tags. Using default \"${CMAKE_GIT_BUILD_ID}\".")
endif()

configure_file(${SRC} ${DST} @ONLY)
