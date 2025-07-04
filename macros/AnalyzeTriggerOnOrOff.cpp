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


static const std::vector<std::string> g_triggersOfInterest = {
    // ---- MBD only ---------------------------------------------------------
    "MBD_NandS_geq_1",
    "MBD_NandS_geq_2",
    "MBD_NandS_geq_2_vtx_lt_10cm",
    "MBD_NandS_geq_2_vtx_lt_30cm",
    "MBD_NandS_geq_2_vtx_lt_60cm",

    // ---- Jet + MBD N&S ≥ 2 -----------------------------------------------
    "Jet_6_GeV_plus_MBD_NS_geq_2",
    "Jet_8_GeV_plus_MBD_NS_geq_2",
    "Jet_10_GeV_plus_MBD_NS_geq_2",
    "Jet_12_GeV_plus_MBD_NS_geq_2",

    // ---- Jet + MBD N&S ≥ 2  with |z_vtx|<10 cm ---------------------------
    "Jet_6_GeV_plus_MBD_NS_geq_2_vtx_lt_10cm",
    "Jet_8_GeV_plus_MBD_NS_geq_2_vtx_lt_10cm",
    "Jet_10_GeV_plus_MBD_NS_geq_2_vtx_lt_10cm",
    "Jet_12_GeV_plus_MBD_NS_geq_2_vtx_lt_10cm",

    // ---- Photon + MBD N&S ≥ 2 --------------------------------------------
    "Photon_2_GeV_plus_MBD_NS_geq_2",
    "Photon_3_GeV_plus_MBD_NS_geq_2",
    "Photon_4_GeV_plus_MBD_NS_geq_2",
    "Photon_5_GeV_plus_MBD_NS_geq_2",

    // ---- Photon + MBD N&S ≥ 2  with |z_vtx|<10 cm ------------------------
    "Photon_3_GeV_plus_MBD_NS_geq_2_vtx_lt_10cm",
    "Photon_4_GeV_plus_MBD_NS_geq_2_vtx_lt_10cm",
    "Photon_5_GeV_plus_MBD_NS_geq_2_vtx_lt_10cm"
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
// 6)  FAST combination analysis
//     – enumerates *only* the combinations that actually occur.
//     – uses a 64‑bit mask instead of std::set keys.
//     – MAX_ORDER == 0 → no limit (all orders); otherwise limit coincidence size.
// --------------------------------------------------------------------------
void analyzeCombinationsFromCSV(const std::string& csvFilePath,
                                unsigned int MAX_ORDER = 3)   // ← tweak here
{
    const std::vector<std::string>& triggers = g_triggersOfInterest;
    const std::size_t N = triggers.size();
    if (N > 64) {
        std::cerr << "[ERROR] Too many triggers (" << N
                  << ") for 64‑bit bitmask. Split the analysis.\n";
        return;
    }

    // -------------------------------------------------
    // 1.  Build runNumber → bitmask of ON triggers
    // -------------------------------------------------
    std::ifstream file(csvFilePath);
    if (!file.is_open()) { std::cerr << "[ERROR] Cannot open " << csvFilePath << "\n"; return; }

    std::string line;
    if (!std::getline(file, line)) { std::cerr << "[ERROR] Empty CSV\n"; return; }

    /* header → column index */
    std::vector<std::string> hdr;
    { std::istringstream iss(line); std::string tok;
      while (std::getline(iss, tok, ',')) hdr.emplace_back(tok); }

    int runCol = -1;
    std::vector<int> trigCol(N, -1);
    for (std::size_t c = 0; c < hdr.size(); ++c) {
        if (hdr[c] == "runNumber") { runCol = static_cast<int>(c); continue; }
        for (std::size_t i = 0; i < N; ++i)
            if (hdr[c] == triggers[i]) trigCol[i] = static_cast<int>(c);
    }
    if (runCol < 0) { std::cerr << "[ERROR] runNumber column missing\n"; return; }

    struct ComboInfo { std::size_t count{}; std::vector<int> runs; };
    std::unordered_map<uint64_t, ComboInfo> combo;   // key = bitmask
    combo.reserve(4096);

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::vector<std::string> tok; tok.reserve(hdr.size());
        std::istringstream iss(line); std::string t;
        while (std::getline(iss, t, ',')) tok.emplace_back(t);

        const int runNo = std::stoi(tok[runCol]);

        uint64_t mask = 0;
        for (std::size_t i = 0; i < N; ++i)
            if (trigCol[i] >= 0 && tok[trigCol[i]] == "ON")
                mask |= (1ULL << i);

        /* iterate over all non‑zero sub‑sets of mask */
        for (uint64_t sub = mask; sub; sub = (sub - 1) & mask) {
#ifdef __cpp_lib_int_pow2     /* C++20 popcount */
            unsigned bits = std::popcount(sub);
#else                           /* GCC/Clang builtin */
            unsigned bits = __builtin_popcountll(sub);
#endif
            if (MAX_ORDER && bits > MAX_ORDER) continue;

            auto &info = combo[sub];
            ++info.count;
            info.runs.push_back(runNo);
        }
    }
    file.close();

    // -------------------------------------------------
    // 2.  Dump statistics (sorted by subset size then lexicographically)
    // -------------------------------------------------
    std::vector<std::pair<uint64_t, ComboInfo>> vec(combo.begin(), combo.end());
    std::sort(vec.begin(), vec.end(), [](auto &a, auto &b) {
#ifdef __cpp_lib_int_pow2
        unsigned sa = std::popcount(a.first), sb = std::popcount(b.first);
#else
        unsigned sa = __builtin_popcountll(a.first), sb = __builtin_popcountll(b.first);
#endif
        return (sa != sb) ? sa < sb : a.first < b.first;
    });

    std::cout << "\nSummary of Trigger Combinations"
              << (MAX_ORDER ? " (order ≤ " + std::to_string(MAX_ORDER) + ")" : "")
              << ":\n";
    for (auto &kv : vec) {
        const uint64_t m = kv.first;
        const ComboInfo &inf = kv.second;

        std::cout << "Combination: ";
        for (std::size_t i = 0; i < N; ++i)
            if (m & (1ULL << i)) std::cout << triggers[i] << " ";
        std::cout << "\nNumber of runs: " << inf.count << "\nRun numbers: ";
        for (std::size_t i = 0; i < inf.runs.size(); ++i) {
            std::cout << inf.runs[i] << " ";
            if ((i + 1) % 10 == 0) std::cout << "\n             ";
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
