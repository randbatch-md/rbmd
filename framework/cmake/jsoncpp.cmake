#安装json cpp方法
#version: jsoncpp-1.9.5.tar.gz

set(JSON_PKG_NAME jsoncpp-1.9.5)
set(JSON_PKG ${JSON_PKG_NAME}.tar.gz)

#json包的路径
set(TOOLS_PATH ${CMAKE_CURRENT_LIST_DIR}/../tools)
set(JSON_PKG_PATH ${TOOLS_PATH}/${JSON_PKG})

#不同的系统安装在不同的目录下面
set(JSON_INSTALL_PATH ${TOOLS_PATH}/${CMAKE_SYSTEM_NAME}/${CMAKE_BUILD_TYPE}/jsoncpp)

#定义安装json宏
macro(setup_jsoncpp)
	#避免重复安装
	if(NOT EXISTS ${JSON_INSTALL_PATH})
		#解压压缩文件
		message("tar xf jsoncpp-1.9.5.tar.gz")
		execute_process(COMMAND ${CMAKE_COMMAND}
		-E tar xf ${JSON_PKG_PATH}
		WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
	
		#cmake配置 cmake -S . -B build -DCMAKE_INSTALL_PREFIX=${JSON_INSTALL_PATH}
		message("cmake -S . -B build jsoncpp-1.9.5")
		set(JSON_SOURCE ${PROJECT_BINARY_DIR}/${JSON_PKG_NAME})
		execute_process(COMMAND ${CMAKE_COMMAND} -S ${JSON_SOURCE} -B ${JSON_SOURCE}/build
		-DCMAKE_INSTALL_PREFIX=${JSON_INSTALL_PATH} -DBUILD_SHARED_LIBS=false)
	
		#编译 cmake --build build
		message("cmake --build build jsoncpp-1.9.5")
		execute_process(COMMAND ${CMAKE_COMMAND} --build ${JSON_SOURCE}/build
		--config ${CMAKE_BUILD_TYPE})
	
		#安装 json
		message("cmake --install jsoncpp-1.9.5")
		execute_process(COMMAND ${CMAKE_COMMAND} --install ${JSON_SOURCE}/build
		--config ${CMAKE_BUILD_TYPE})
	endif()
endmacro()