function [Air1] = importfile_diffspecies(workbookFile, sheetName, dataLines)
%IMPORTFILE Import data from a spreadsheet
%  [AIR1, VARNAME2, VARNAME3, VARNAME4, VARNAME5, VARNAME6, VARNAME7,
%  VARNAME8, VARNAME9, VARNAME10, VARNAME11, VARNAME12, VARNAME13,
%  VARNAME14, VARNAME15, VARNAME16, VARNAME17, VARNAME18, VARNAME19,
%  VARNAME20, VARNAME21, VARNAME22, VARNAME23, VARNAME24, VARNAME25] =
%  IMPORTFILE(FILE) reads data from the first worksheet in the Microsoft
%  Excel spreadsheet file named FILE.  Returns the data as column
%  vectors.
%
%  [AIR1, VARNAME2, VARNAME3, VARNAME4, VARNAME5, VARNAME6, VARNAME7,
%  VARNAME8, VARNAME9, VARNAME10, VARNAME11, VARNAME12, VARNAME13,
%  VARNAME14, VARNAME15, VARNAME16, VARNAME17, VARNAME18, VARNAME19,
%  VARNAME20, VARNAME21, VARNAME22, VARNAME23, VARNAME24, VARNAME25] =
%  IMPORTFILE(FILE, SHEET) reads from the specified worksheet.
%
%  [AIR1, VARNAME2, VARNAME3, VARNAME4, VARNAME5, VARNAME6, VARNAME7,
%  VARNAME8, VARNAME9, VARNAME10, VARNAME11, VARNAME12, VARNAME13,
%  VARNAME14, VARNAME15, VARNAME16, VARNAME17, VARNAME18, VARNAME19,
%  VARNAME20, VARNAME21, VARNAME22, VARNAME23, VARNAME24, VARNAME25] =
%  IMPORTFILE(FILE, SHEET, DATALINES) reads from the specified worksheet
%  for the specified row interval(s). Specify DATALINES as a positive
%  scalar integer or a N-by-2 array of positive scalar integers for
%  dis-contiguous row intervals.
%
%  Example:
%  [Air1, VarName2, VarName3, VarName4, VarName5, VarName6, VarName7, VarName8, VarName9, VarName10, VarName11, VarName12, VarName13, VarName14, VarName15, VarName16, VarName17, VarName18, VarName19, VarName20, VarName21, VarName22, VarName23, VarName24, VarName25] = importfile("C:\Users\vdela\OneDrive\Desktop\Cantera_VDL\Djk.xlsx", "Sheet1", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 29-Aug-2024 14:47:59

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [1, Inf];
end

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 25);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = dataLines(1, :);

% Specify column names and types
opts.VariableNames = ["Air1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", "VarName20", "VarName21", "VarName22", "VarName23", "VarName24", "VarName25"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "Air1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Air1", "EmptyFieldRule", "auto");

% Import the data
tbl = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = dataLines(idx, :);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    tbl = [tbl; tb]; %#ok<AGROW>
end

%% Convert to output type
Air1 = tbl.Air1;
VarName2 = tbl.VarName2;
VarName3 = tbl.VarName3;
VarName4 = tbl.VarName4;
VarName5 = tbl.VarName5;
VarName6 = tbl.VarName6;
VarName7 = tbl.VarName7;
VarName8 = tbl.VarName8;
VarName9 = tbl.VarName9;
VarName10 = tbl.VarName10;
VarName11 = tbl.VarName11;
VarName12 = tbl.VarName12;
VarName13 = tbl.VarName13;
VarName14 = tbl.VarName14;
VarName15 = tbl.VarName15;
VarName16 = tbl.VarName16;
VarName17 = tbl.VarName17;
VarName18 = tbl.VarName18;
VarName19 = tbl.VarName19;
VarName20 = tbl.VarName20;
VarName21 = tbl.VarName21;
VarName22 = tbl.VarName22;
VarName23 = tbl.VarName23;
VarName24 = tbl.VarName24;
VarName25 = tbl.VarName25;
end