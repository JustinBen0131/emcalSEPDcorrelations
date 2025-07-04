#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdio>  // for std::remove()

// ROOT
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

// ---------------------------------------------------------------
// 1) The DB -> Folder map (LHS -> RHS)
//    And the CSV columns in the order of the RHS strings
// ---------------------------------------------------------------
static std::map<std::string,std::string> g_dbNameToFolderName = {
    {"MBD N&S >= 1",              "MBD_NandS_geq_1"},
    {"MBD N&S >= 2",              "MBD_NandS_geq_2"},
    {"MBD N&S >= 2, vtx < 10 cm", "MBD_NandS_geq_2_vtx_lt_10cm"},
    {"MBD N&S >= 2, vtx < 30 cm", "MBD_NandS_geq_2_vtx_lt_30cm"},
    {"MBD N&S >= 2, vtx < 60 cm", "MBD_NandS_geq_2_vtx_lt_60cm"},
    {"Jet 6 GeV + MBD NS >=2",    "Jet_6_GeV_plus_MBD_NS_geq_2"},
    {"Jet 8 GeV + MBD NS >= 2",   "Jet_8_GeV_plus_MBD_NS_geq_2"},
    {"Jet 10 GeV + MBD NS >= 2",  "Jet_10_GeV_plus_MBD_NS_geq_2"},
    {"Jet 12 GeV + MBD NS >= 2",  "Jet_12_GeV_plus_MBD_NS_geq_2"},
    {"Photon 2 GeV+ MBD NS >= 2", "Photon_2_GeV_plus_MBD_NS_geq_2"},
    {"Photon 3 GeV + MBD NS >= 2","Photon_3_GeV_plus_MBD_NS_geq_2"},
    {"Photon 4 GeV + MBD NS >= 2","Photon_4_GeV_plus_MBD_NS_geq_2"},
    {"Photon 5 GeV + MBD NS >= 2","Photon_5_GeV_plus_MBD_NS_geq_2"},
    {"Jet 6 GeV + MBD NS >=2, vtx < 10 cm", "Jet_6_GeV_plus_MBD_NS_geq_2_vtx_lt_10cm"},
    {"Jet 8 GeV + MBD NS >=2, vtx < 10 cm", "Jet_8_GeV_plus_MBD_NS_geq_2_vtx_lt_10cm"},
    {"Jet 10 GeV + MBD NS >=2, vtx < 10 cm", "Jet_10_GeV_plus_MBD_NS_geq_2_vtx_lt_10cm"},
    {"Jet 12 GeV + MBD NS >=2, vtx < 10 cm", "Jet_12_GeV_plus_MBD_NS_geq_2_vtx_lt_10cm"},
    {"Photon 3 GeV + MBD NS >= 2, vtx < 10 cm","Photon_3_GeV_plus_MBD_NS_geq_2_vtx_lt_10cm"},
    {"Photon 4 GeV + MBD NS >= 2, vtx < 10 cm","Photon_4_GeV_plus_MBD_NS_geq_2_vtx_lt_10cm"},
    {"Photon 5 GeV + MBD NS >= 2, vtx < 10 cm","Photon_5_GeV_plus_MBD_NS_geq_2_vtx_lt_10cm"}
};

// The exact order we want in the CSV
static const std::vector<std::string> g_triggersOfInterest = {
    "MBD_NandS_geq_1",
    "Jet_6_GeV_plus_MBD_NS_geq_1",
    "Jet_8_GeV_plus_MBD_NS_geq_1",
    "Jet_10_GeV_plus_MBD_NS_geq_1",
    "Jet_12_GeV_plus_MBD_NS_geq_1",
    "Photon_2_GeV_plus_MBD_NS_geq_1",
    "Photon_3_GeV_plus_MBD_NS_geq_1",
    "Photon_4_GeV_plus_MBD_NS_geq_1",
    "Photon_5_GeV_plus_MBD_NS_geq_1"
};

// --------------------------------------------------------------------------
// 2) LOAD RUN NUMBERS FROM A TEXT FILE
// --------------------------------------------------------------------------
std::vector<int> loadRunNumbersFromFile(const std::string& filename)
{
    std::vector<int> runNumbers;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "[ERROR] Cannot open run list file: " << filename << std::endl;
        return runNumbers;
    }
    std::string line;
    while (std::getline(file, line)) {
        try {
            int runNumber = std::stoi(line);
            runNumbers.push_back(runNumber);
        }
        catch (const std::invalid_argument&) {
            std::cerr << "[WARNING] Invalid run number in file: " << line << std::endl;
        }
    }
    file.close();
    return runNumbers;
}

// --------------------------------------------------------------------------
// 3) For a given runNumber, query gl1_scalers JOIN gl1_triggernames
//    Collect which DB triggers had scaled > 0
//    Then map them to the "folder name" via g_dbNameToFolderName
//    Return: folderName => ON(true)/OFF(false)
// --------------------------------------------------------------------------
std::unordered_map<std::string,bool> getFolderOnOffForRun(int runNumber)
{
    std::unordered_map<std::string,bool> folderOnOff;

    // Connect to DB
    TSQLServer* db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq", "phnxro", "");
    if (!db || db->IsZombie()) {
        std::cerr << "[ERROR] DB connection failed for run " << runNumber << std::endl;
        if (db) delete db;
        return folderOnOff;
    }

    // Build query: get (t.triggername, s.scaled)
    char query[1024];
    snprintf(query, sizeof(query),
        "SELECT t.triggername, s.scaled "
        " FROM gl1_scalers s "
        " JOIN gl1_triggernames t ON s.index = t.index "
        "    AND s.runnumber BETWEEN t.runnumber AND t.runnumber_last "
        " WHERE s.runnumber = %d;",
        runNumber);

    TSQLResult* res = db->Query(query);
    if (!res) {
        std::cerr << "[ERROR] Query failed for run " << runNumber << std::endl;
        delete db;
        return folderOnOff;
    }

    // For each row => check if scaled>0 => ON => store in folderOnOff
    while (TSQLRow* row = res->Next()) {
        const char* dbTriggerName = row->GetField(0);
        const char* scaledStr     = row->GetField(1);
        if (!dbTriggerName || !scaledStr) {
            delete row;
            continue;
        }
        std::string trigName(dbTriggerName);
        int scaledCounts = std::stoi(scaledStr);
        bool isOn = (scaledCounts > 0);

        // Now see if this DB trigger name is in g_dbNameToFolderName
        auto mapIt = g_dbNameToFolderName.find(trigName);
        if (mapIt != g_dbNameToFolderName.end()) {
            // If found => the folder name is mapIt->second
            std::string folderName = mapIt->second;
            folderOnOff[folderName] = isOn;
        }
        delete row;
    }

    delete res;
    delete db;
    return folderOnOff;
}

// --------------------------------------------------------------------------
// 4) For each run, query DB => build a map folderName => ON/OFF
//    Then write out CSV row: runNumber, col1..col9 => ON/OFF
//    BEFORE writing, remove any old CSV file with the same name
// --------------------------------------------------------------------------
void analyzeTriggersAndWriteCSV_DB(const std::vector<int>& runNumbers,
                                   const std::string& outputCSV)
{
    // Delete any existing CSV with this name
    std::remove(outputCSV.c_str());

    std::ofstream csvFile(outputCSV);
    if (!csvFile.is_open()) {
        std::cerr << "[ERROR] Cannot open output CSV: " << outputCSV << std::endl;
        return;
    }

    // Write header
    csvFile << "runNumber";
    for (const auto& folderName : g_triggersOfInterest) {
        csvFile << "," << folderName;
    }
    csvFile << "\n";

    // For each run, do DB query => fill row
    for (int runNum : runNumbers) {
        // Query DB => get folder->bool
        std::unordered_map<std::string,bool> folderOnOff = getFolderOnOffForRun(runNum);

        // Start row with runNumber
        csvFile << runNum;

        // For each folder in our order => check if ON or OFF
        for (const auto& folderName : g_triggersOfInterest) {
            bool isOn = false;
            auto it = folderOnOff.find(folderName);
            if (it != folderOnOff.end()) {
                isOn = it->second;
            }
            csvFile << (isOn ? ",ON" : ",OFF");
        }
        csvFile << "\n";
    }

    csvFile.close();
    std::cout << "[INFO] CSV written: " << outputCSV << std::endl;
}

// --------------------------------------------------------------------------
// 5) Summarize Trigger Counts from CSV
// --------------------------------------------------------------------------
void summarizeTriggerCountsFromCSV(const std::string& csvFilePath,
                                   const std::vector<std::string>& triggersOfInterest)
{
    std::ifstream inFile(csvFilePath);
    if (!inFile.is_open()) {
        std::cerr << "[ERROR] Cannot open CSV: " << csvFilePath << std::endl;
        return;
    }

    // Skip header
    std::string header;
    if (!std::getline(inFile, header)) {
        std::cerr << "[ERROR] CSV is empty: " << csvFilePath << std::endl;
        return;
    }

    // ON/OFF counters
    std::unordered_map<std::string,int> onCount, offCount;
    for (auto& trig : triggersOfInterest) {
        onCount[trig] = 0;
        offCount[trig] = 0;
    }

    // Read lines
    std::string line;
    while (std::getline(inFile, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;
        while (std::getline(iss, token, ',')) {
            // trim
            token.erase(0, token.find_first_not_of(" \t\r\n"));
            token.erase(token.find_last_not_of(" \t\r\n") + 1);
            tokens.push_back(token);
        }
        // Expect: runNumber + #triggers
        if (tokens.size() != (1 + triggersOfInterest.size())) continue;

        // tokens[0] = runNumber, skip that
        for (size_t i=0; i<triggersOfInterest.size(); i++) {
            std::string val = tokens[i+1]; // offset by 1
            if (val=="ON")  onCount[triggersOfInterest[i]]++;
            else            offCount[triggersOfInterest[i]]++;
        }
    }
    inFile.close();

    // Print summary
    std::cout << "\nSummary of Trigger Counts:\n";
    for (auto& trig : triggersOfInterest) {
        std::cout << "  " << trig << ":\n"
                  << "    ON  = " << onCount[trig] << "\n"
                  << "    OFF = " << offCount[trig] << "\n"
                  << "----------------------------------------\n";
    }
}

// --------------------------------------------------------------------------
// 6) Optional: Analyze Combinations from the CSV
// --------------------------------------------------------------------------
void analyzeCombinationsFromCSV(const std::string& csvFilePath)
{
    // We'll do combos of all 9 triggers in g_triggersOfInterest
    std::vector<std::string> triggers = g_triggersOfInterest;

    // Read CSV
    std::ifstream file(csvFilePath);
    if (!file.is_open()) {
        std::cerr << "[ERROR] Cannot open CSV for combos: " << csvFilePath << std::endl;
        return;
    }
    std::string headerLine;
    if (!std::getline(file, headerLine)) {
        std::cerr << "[ERROR] CSV is empty: " << csvFilePath << std::endl;
        return;
    }

    // Parse header => find columns
    std::vector<std::string> headers;
    {
        std::istringstream hdr(headerLine);
        std::string h;
        while (std::getline(hdr, h, ',')) {
            h.erase(0, h.find_first_not_of(" \t\r\n"));
            h.erase(h.find_last_not_of(" \t\r\n") + 1);
            headers.push_back(h);
        }
    }

    int runNumberCol = -1;
    std::map<std::string,int> trigToCol;
    for (size_t c=0; c<headers.size(); c++) {
        if (headers[c]=="runNumber") {
            runNumberCol = c;
        }
        else {
            // If it's in triggers, store it
            for (auto& t: triggers) {
                if (headers[c]==t) {
                    trigToCol[t] = (int)c;
                }
            }
        }
    }
    if (runNumberCol<0) {
        std::cerr << "[ERROR] 'runNumber' column not found in CSV." << std::endl;
        return;
    }

    // We'll store run => trig => ON/OFF
    std::map<int,std::map<std::string,std::string>> runTrigStatus;
    std::vector<int> allRuns;

    // Read each line
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string t;
        while (std::getline(iss, t, ',')) {
            t.erase(0, t.find_first_not_of(" \t\r\n"));
            t.erase(t.find_last_not_of(" \t\r\n") + 1);
            tokens.push_back(t);
        }
        if (tokens.size()!=headers.size()) continue;
        int runNo = std::stoi(tokens[runNumberCol]);
        allRuns.push_back(runNo);

        for (auto& trName : triggers) {
            auto it = trigToCol.find(trName);
            if (it==trigToCol.end()) {
                // column not found => OFF
                runTrigStatus[runNo][trName] = "OFF";
            } else {
                runTrigStatus[runNo][trName] = tokens[it->second];
            }
        }
    }
    file.close();

    // Generate combos
    int n = (int)triggers.size();
    int total = (1<<n);
    std::map<std::set<std::string>,std::vector<int>> comboToRuns;

    for (int mask=1; mask<total; mask++) {
        std::set<std::string> combo;
        for (int i=0; i<n; i++) {
            if (mask&(1<<i)) {
                combo.insert(triggers[i]);
            }
        }
        // check each run => if all triggers in combo == ON
        for (int rn : allRuns) {
            bool allOn = true;
            for (auto& trig : combo) {
                if (runTrigStatus[rn][trig]!="ON") {
                    allOn=false;
                    break;
                }
            }
            if (allOn) {
                comboToRuns[combo].push_back(rn);
            }
        }
    }

    // Sort combos by size
    std::vector<std::pair<std::set<std::string>,std::vector<int>>> sorted;
    for (auto& kv : comboToRuns) {
        sorted.push_back(kv);
    }
    std::sort(sorted.begin(), sorted.end(),
              [](auto& a, auto& b){
                  if (a.first.size()!=b.first.size()) {
                      return a.first.size()<b.first.size();
                  }
                  return a.first < b.first;
              });

    // Print
    std::cout << "\nSummary of Trigger Combinations:\n";
    for (auto& kv : sorted) {
        const auto& combo = kv.first;
        const auto& runs  = kv.second;
        std::cout << "Combination: ";
        for (auto& c : combo) {
            std::cout << c << " ";
        }
        std::cout << "\nNumber of runs: " << runs.size() << "\nRun numbers: ";
        for (size_t i=0; i<runs.size(); i++) {
            std::cout << runs[i] << " ";
            if ((i+1)%10==0) std::cout << "\n             ";
        }
        std::cout << "\n-------------------------------------\n";
    }
}

// --------------------------------------------------------------------------
// 7) MAIN "AnalyzeTriggerOnOrOff" style function
// --------------------------------------------------------------------------
void AnalyzeTriggerOnOrOff()
{
    // 1) Load from a single text file
    std::string inputFile = "/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/goodRunList_sEPD_run24auau.txt";
    std::vector<int> runNumbers = loadRunNumbersFromFile(inputFile);
    if (runNumbers.empty()) {
        std::cerr << "[ERROR] No runs found in: " << inputFile << std::endl;
        return;
    }
    std::cout << "[INFO] Loaded " << runNumbers.size() << " runs.\n";

    // 2) Output CSV path: "triggerAnalysisCombined.csv"
    std::string csvOut = "/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/triggerOnOFF.csv";

    // 3) Build the CSV (removing old file if present)
    analyzeTriggersAndWriteCSV_DB(runNumbers, csvOut);

    // 4) Summaries
    summarizeTriggerCountsFromCSV(csvOut, g_triggersOfInterest);

    // 5) Combinations
    analyzeCombinationsFromCSV(csvOut);
}
