function udisp(str,guiRunning,guiHandle)

if guiRunning
 guiHandle.gui.eto.myEstModelPopup.append(str);
else
 disp(str);
end