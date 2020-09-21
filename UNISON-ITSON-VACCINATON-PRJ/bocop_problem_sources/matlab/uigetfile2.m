%Copyright (c) 2009, Chris Cannell
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without 
%modification, are permitted provided that the following conditions are 
%met:
%
%    * Redistributions of source code must retain the above copyright 
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright 
%      notice, this list of conditions and the following disclaimer in 
%      the documentation and/or other materials provided with the distribution
%    * Neither the name of the Virginia Tech nor the names 
%      of its contributors may be used to endorse or promote products derived 
%      from this software without specific prior written permission.
%      
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%POSSIBILITY OF SUCH DAMAGE.




function [filename, pathname, filterindex] = uigetfile2(varargin)

%UIGETFILE2 Standard open file dialog box which remembers last opened folder
%   UIGETFILE2 is a wrapper for Matlab's UIGETFILE function which adds the
%   ability to remember the last folder opened.  UIGETFILE2 stores
%   information about the last folder opened in a mat file which it looks
%   for when called.
%
%   UIGETFILE2 can only remember the folder used if the current directory
%   is writable so that a mat file can be stored.  Only successful file
%   selections update the folder remembered.  If the user cancels the file
%   dialog box then the remembered path is left the same.
%
%   Usage is the same as UIGETFILE.
%
%
%   See also UIGETFILE, UIPUTFILE, UIGETDIR.

%   Written by Chris J Cannell and Aditya Gadre
%   Contact ccannell@mindspring.com for questions or comments.
%   12/05/2005

% name of mat file to save last used directory information
lastDirMat = 'lastUsedDir.mat';

% save the present working directory
savePath = pwd;
% set default dialog open directory to the present working directory
lastDir = savePath;
% load last data directory
if exist(lastDirMat, 'file') ~= 0
    % lastDirMat mat file exists, load it
    load('-mat', lastDirMat)
    % check if lastDataDir variable exists and contains a valid path
    if (exist('lastUsedDir', 'var') == 1) && ...
            (exist(lastUsedDir, 'dir') == 7)
        % set default dialog open directory
        lastDir = lastUsedDir;
    end
end

% load folder to open dialog box in
cd(lastDir);
% call uigetfile with arguments passed from uigetfile2 function
[filename, pathname, filterindex] = uigetfile(varargin{:});
% change path back to original working folder
cd(savePath);

% if the user did not cancel the file dialog then update lastDirMat mat
% file with the folder used
if ~isequal(filename,0) && ~isequal(pathname,0)
    try
        % save last folder used to lastDirMat mat file
        lastUsedDir = pathname;
        save(lastDirMat, 'lastUsedDir');
    catch
        % error saving lastDirMat mat file, display warning, the folder
        % will not be remembered
        disp(['Warning: Could not save file ''', lastDirMat, '''']);
    end
end
