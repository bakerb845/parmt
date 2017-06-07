add_sources(
angleMT.c
beta2u.c
cart2sph.c
CMT2m0.c
CMT2mw.c
convertMT.c
eulerUtil_rotmat.c
eulerUtil_rotmatGen.c
gamma2v.c
h2theta.c
lam2phizeta.c
lamsort.c
lune2lam.c
lune2rect.c
m02hdur.c
m02mw.c
mw2m0.c
Mvec2Mmat.c
normal2strdip.c
normMat.c
normMT.c
rect2lune.c
tape2015Eqn7.c
theta2h.c
tt2cmt.c
u2beta.c
v2gamma.c
xyz2tp.c
)
set (SRCS_COMPEARTH ${SRCS} PARENT_SCOPE)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/compearth_config.h)
set(headers compearth.h compearth_config.h)
file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
install(FILES ${headers} DESTINATION ${_relPath} COMPONENT headers)
