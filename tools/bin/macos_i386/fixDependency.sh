#!/bin/sh
echo "Adjusting path of libtwtool in libtwh5"
echo "======================================"
CURRENT_TOOL_PATH=`otool -l libtwh5.dylib | awk '/libtwtool.dylib/{print $2}'`
echo "current tool path: " $CURRENT_TOOL_PATH
NEW_TOOL_PATH="`pwd`/libtwtool.dylib"
echo "    new tool path: " $NEW_TOOL_PATH
install_name_tool -change $CURRENT_TOOL_PATH $NEW_TOOL_PATH libtwh5.dylib
