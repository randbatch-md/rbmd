#安装cxxopts方法
#version: cxxopts-3.1.1.tar.gz

set(CXXOPTS_PKG_NAME cxxopts-3.1.1)
set(CXXOPTS_PKG ${CXXOPTS_PKG_NAME}.tar.gz)

#cxxopts包的路径
set(TOOLS_PATH ${CMAKE_CURRENT_LIST_DIR}/../tools)
set(CXXOPTS_PKG_PATH ${TOOLS_PATH}/${CXXOPTS_PKG})

#不同的系统安装在不同的目录下面
set(CXXOPTS_INSTALL_PATH ${TOOLS_PATH}/${CMAKE_SYSTEM_NAME}/${CMAKE_BUILD_TYPE}/cxxopts)

#定义安装cxxopts宏
macro(setup_cxxopts)
	#避免重复安装
	if(NOT EXISTS ${CXXOPTS_INSTALL_PATH})
		#解压压缩文件
		message("tar xf cxxopts-3.1.1.tar.gz")
		execute_process(COMMAND ${CMAKE_COMMAND}
		-E tar xf ${CXXOPTS_PKG_PATH}
		WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
	
		#cmake配置 cmake -S . -B build -DCMAKE_INSTALL_PREFIX=${CXXOPTS_INSTALL_PATH}
		message("cmake -S . -B build cxxopts-3.1.1")
		set(CXXOPTS_SOURCE ${PROJECT_BINARY_DIR}/${CXXOPTS_PKG_NAME})
		execute_process(COMMAND ${CMAKE_COMMAND} -S ${CXXOPTS_SOURCE} -B ${CXXOPTS_SOURCE}/build
		-DCMAKE_INSTALL_PREFIX=${CXXOPTS_INSTALL_PATH} -DBUILD_SHARED_LIBS=false)
	
		#编译 cmake --build build
		message("cmake --build build cxxopts-3.1.1")
		execute_process(COMMAND ${CMAKE_COMMAND} --build ${CXXOPTS_SOURCE}/build
		--config ${CMAKE_BUILD_TYPE})
	
		#安装 cxxopts
		message("cmake --install cxxopts-3.1.1")
		execute_process(COMMAND ${CMAKE_COMMAND} --install ${CXXOPTS_SOURCE}/build
		--config ${CMAKE_BUILD_TYPE})
	endif()
endmacro()