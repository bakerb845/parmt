configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/prepmt_config.h)
set(headers prepmt_commands.h prepmt_event.h prepmt.h prepmt_hudson96.h prepmt_prepData.h
prepmt_dataArchive.h  prepmt_greens.h  prepmt_hpulse96.h  prepmt_pickFile.h  prepmt_struct.h)
file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
install(FILES ${headers} DESTINATION ${_relPath} COMPONENT headers)
