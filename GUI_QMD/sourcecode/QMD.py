# -*- coding: utf-8 -*-
from QMD_gui import Ui_QMD
import sys
from PyQt5.QtGui import QIntValidator,QDoubleValidator,QDesktopServices
from PyQt5.QtCore import QDir,QUrl
from PyQt5.QtWidgets import QAction, QFileDialog, QMainWindow,QApplication
import os
import re
import pandas as pd
import platform
import multiprocessing
import qmd_dataPreprocessing_exe
import qmd_optimization_exe
import load_save_project
import datetime


# path = getattr(sys, '_MEIPASS', os.getcwd())   
# os.chdir(path)

s = platform.uname()
os_p = s[0]

cpu_count = multiprocessing.cpu_count()


class query_window(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.ui = Ui_QMD()
        self.ui.setupUi(self)
        self.ui.lineEdit_6.setValidator(QIntValidator())
        self.ui.lineEdit_9.setValidator(QDoubleValidator())
        self.ui.lineEdit_4.setValidator(QIntValidator(1, 1024))
        self.ui.pushButton.clicked.connect(self.reset_all_input)
        self.ui.pushButton_2.clicked.connect(self.run_qmd)
        self.ui.pushButton_3.clicked.connect(self.getFiles)
        self.ui.comboBox_4.setCurrentIndex(1)

        self.LoadProject = QAction("Load Project")
        self.LoadProject.setShortcut('Ctrl+o')
        self.LoadProject.triggered.connect(self.Load_Project)
        self.ui.menu.addAction(self.LoadProject)

        self.openmanual = QAction("Open Manual")
        self.openmanual.setShortcut('Ctrl+m')
        self.openmanual.triggered.connect(self.open_manual)
        self.ui.menu.addAction(self.openmanual)

        self.aboudQMD = QAction("About QMD")
        self.aboudQMD.setShortcut('Ctrl+a')
        self.aboudQMD.triggered.connect(self.aboud_QMD)
        self.ui.menu.addAction(self.aboudQMD)

        self.exit = QAction("Exit Application")
        self.exit.setShortcut('Ctrl+q')
        self.exit.triggered.connect(self.exit_app)
        self.ui.menu.addAction(self.exit)

    def Load_Project(self):
        # not done
        dig = QFileDialog()
        dig.setNameFilters(["qmd project file(*.qmd)"])
        dig.setFileMode(QFileDialog.ExistingFile)
        dig.setFilter(QDir.Files)

        if dig.exec_():
            filenames = dig.selectedFiles()
            if os_p == 'Windows':
                tmpdatafilename = filenames[0].replace('/', '\\')
            else:
                tmpdatafilename = filenames[0]
            try:
                project_dict = load_save_project.load_dict(tmpdatafilename)
                self.ui.comboBox_3.setCurrentIndex(project_dict.get('fdr'))
                self.ui.comboBox_4.setCurrentIndex(
                    project_dict.get('confidenceI'))
                self.ui.comboBox_5.setCurrentIndex(project_dict.get('plot'))
                self.ui.lineEdit_2.setText(project_dict.get('projectname'))
                self.ui.lineEdit_3.setText(project_dict.get('datafile'))
                self.ui.lineEdit_4.setText(project_dict.get(
                    'minimum_taxa_detection_samples'))
                self.ui.lineEdit_5.setText(project_dict.get('grouplabel'))
                self.ui.lineEdit_6.setText(project_dict.get('permu'))
                self.ui.lineEdit_7.setText(project_dict.get('treatlabel'))
                self.ui.lineEdit_8.setText(project_dict.get('controllabel'))
                self.ui.lineEdit_9.setText(
                    project_dict.get('minimum_taxa_median_abundance'))
                self.ui.textEdit.setText(project_dict.get('loggs'))
            except:
                self.ui.textEdit.setText('Project recovery failed.')

    def getFiles(self):
        dig = QFileDialog()
        dig.setNameFilters(["csv data file(*.csv)"])
        dig.setFileMode(QFileDialog.ExistingFile)
        dig.setFilter(QDir.Files)
        if dig.exec_():
            filenames = dig.selectedFiles()
            if os_p == 'Windows':
                tmpdatafilename = filenames[0].replace('/', '\\')
            else:
                tmpdatafilename = filenames[0]
            self.ui.lineEdit_3.setText(tmpdatafilename)

    def aboud_QMD(self):
        file_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "source/About_QMD.html"))
        local_url = QUrl.fromLocalFile(file_path)
        QDesktopServices.openUrl(local_url)


    def open_manual(self):
        file_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "source/manual.html"))
        local_url = QUrl.fromLocalFile(file_path)
        QDesktopServices.openUrl(local_url)

    def exit_app(self):
        self.close()

    def menuquit(self):
        self.close()

    def reset_all_input(self):
        self.ui.comboBox_3.setCurrentIndex(0)
        self.ui.comboBox_4.setCurrentIndex(1)
        self.ui.comboBox_5.setCurrentIndex(0)
        self.ui.lineEdit_2.setText('your project name')
        self.ui.lineEdit_3.setText('D:\\data\\abundance_to_be_analyzed.csv')
        self.ui.lineEdit_4.setText('5')
        self.ui.lineEdit_5.setText('condition')
        self.ui.lineEdit_6.setText('500')
        self.ui.lineEdit_7.setText('treat')
        self.ui.lineEdit_8.setText('control')
        self.ui.lineEdit_9.setText('0')
        self.ui.textEdit.setText('Logs..')

    def run_qmd(self):
        if os_p == 'Windows':
            datafilename = self.ui.lineEdit_3.text().replace('/', '\\')
        else:
            datafilename = self.ui.lineEdit_3.text()

        projectname = re.sub("[^a-zA-Z0-9]", "_", self.ui.lineEdit_2.text())
        sharp_str = '######################\r######################\r\r'
        currenttime = str(datetime.datetime.now()) + '\r'
        setting_str = sharp_str+currenttime+'Setting:\r'
        setting_str = setting_str+'\tYour project name: '+self.ui.lineEdit_2.text()+'\r' +\
            '\tPermutation Loops: '+self.ui.lineEdit_6.text()+'\r' +\
            '\tConfidence Interval : '+self.ui.comboBox_4.currentText()+'\r' + \
            '\tFDR adjustment : ' + self.ui.comboBox_3.currentText() + '\r' + \
            '\tPlot : ' + self.ui.comboBox_5.currentText() + '\r' +\
            '\tMinimum taxa detection samples in each group: ' + self.ui.lineEdit_4.text() + '\r' + \
            '\tMinimum taxa median abundance in control group: ' + self.ui.lineEdit_9.text() + '\r' + \
            '\tColumn name of group indicator: ' + self.ui.lineEdit_5.text() + '\r' + \
            '\tLabel of control group: ' + self.ui.lineEdit_8.text() + '\r' + \
            '\tLabel of treatment group: ' + self.ui.lineEdit_7.text() + '\r' + \
            '\tData file: ' + datafilename + '\r'
        QApplication.processEvents()
        self.ui.textEdit.setText(setting_str)
        QApplication.processEvents()
        isExists = os.path.exists(self.ui.lineEdit_3.text())
        if not isExists:
            setting_str = setting_str + '\tThe data file not existed. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if isExists:
            try:
                relativedata = pd.read_csv(
                    datafilename, header=0, index_col=False)
            except:
                setting_str = setting_str + \
                    '\tCan not read into the data file. Please check you data file. Make sure it follows the data template. Quit the analysis.' + '\r'
                self.ui.textEdit.setText(setting_str)
                return None

        # check the parameter

        try:
            group = relativedata[self.ui.lineEdit_5.text()]
        except:
            setting_str = setting_str + '\tCan not find the group indicator column in the data file. Please check you data file. Make sure it follows the data template. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        controlnum = len(
            relativedata[relativedata[self.ui.lineEdit_5.text()] == self.ui.lineEdit_8.text()])
        treatnum = len(
            relativedata[relativedata[self.ui.lineEdit_5.text()] == self.ui.lineEdit_7.text()])
        if controlnum == 0:
            setting_str = setting_str + '\tData file do not contain the control group data. Please check you data file. Make sure it follows the data template. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if treatnum == 0:
            setting_str = setting_str + '\tData file do not contain the treat group data. Please check you data file. Make sure it follows the data template. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if controlnum < 5:
            setting_str = setting_str + \
                '\tThe sample size of control group data is to small. It is preferably larger than 5. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if treatnum < 5:
            setting_str = setting_str + \
                '\tThe sample size of control group data is to small. It is preferably larger than 5. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        minimum_taxa_detection_num = int(self.ui.lineEdit_4.text())
        if minimum_taxa_detection_num < 2:
            setting_str = setting_str + \
                '\tThe minimum taxa detection num should be larger than 2. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if controlnum < minimum_taxa_detection_num:
            setting_str = setting_str + \
                '\tThe minimum taxa detection num should be smaller than the sample size of control group ('+str(
                    controlnum)+'). Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if treatnum < minimum_taxa_detection_num:
            setting_str = setting_str + \
                '\tThe minimum taxa detection num should be smaller than the sample size of treat group ('+str(
                    treatnum)+'). Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None

        tmptaxalist = list(relativedata.columns)
        taxalist = []
        all_float_dtype=True
        countneg=0
        for g in tmptaxalist:
            if g != self.ui.lineEdit_5.text():
                taxalist.append(g)
                if not pd.api.types.is_numeric_dtype(relativedata[g]):
                    all_float_dtype=False
                else:
                    countneg=countneg+sum(relativedata[g]<0)
        if countneg>0:
            setting_str = setting_str + \
                '\tThe data file contains negative data. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if not all_float_dtype:
            setting_str = setting_str + \
                '\tSome columns of the data file are not abundance data. Please check you data file. Make sure it follows the data template. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if len(taxalist) == 0:
            setting_str = setting_str + \
                '\tNo taxa are found in the data file. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if len(taxalist) < 3:
            setting_str = setting_str + '\tOnly ' + \
                str(len(taxalist)) + \
                ' taxa are found in the data file. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None

        setting_str = setting_str + '\tThe data file contains ' + \
            str(len(taxalist))+' taxa' + '\r'
        self.ui.textEdit.setText(setting_str)
        setting_str = setting_str + '\tThe '+self.ui.lineEdit_8.text()+' group contains ' + \
            str(controlnum) + ' samples' + '\r'
        self.ui.textEdit.setText(setting_str)
        setting_str = setting_str + '\tThe ' + self.ui.lineEdit_7.text() + ' group contains ' + \
            str(treatnum) + ' samples' + '\r'
        self.ui.textEdit.setText(setting_str)
        QApplication.processEvents()
        # print(self.ui.lineEdit.text())
        # print(self.ui.comboBox_2.currentText())

        if float(self.ui.lineEdit_9.text()) < 0:
            setting_str = setting_str + \
                '\tThe parameter minimum taxa median abundance in control group should be non-negative. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        currentpath = os.getcwd()+'/Project'
        isExists = os.path.exists(currentpath)
        if not isExists:
            try:
                os.makedirs(currentpath)
            except:
                setting_str = setting_str + \
                    '\tCreating project analysis result storage folder ("' + \
                    currentpath + '") failed. Quit the analysis.' + '\r'
                self.ui.textEdit.setText(setting_str)
                return None
        if os_p == 'Windows':
            fileplace = currentpath.replace(
                '/', '\\').rstrip("\\")+"\\"+projectname
        else:
            fileplace = currentpath.rstrip("/") + "/" + projectname
        isExists = os.path.exists(fileplace)
        if isExists:
            setting_str = setting_str+'\tProject name existed. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if not isExists:
            try:
                os.makedirs(fileplace)
                setting_str = setting_str + \
                    '\tProject analysis result are stored in "' + fileplace + '"\r'
                self.ui.textEdit.setText(setting_str)
                QApplication.processEvents()
            except:
                setting_str = setting_str + \
                    '\tCreating project analysis result storage folder ("' + \
                    fileplace + '") failed. Quit the analysis.' + '\r'
                self.ui.textEdit.setText(setting_str)
                return None
            try:
                if os_p == 'Windows':
                    data_preprocessing_place = fileplace+"\\datapreprocessing"
                else:
                    data_preprocessing_place = fileplace + "/datapreprocessing"
                os.makedirs(data_preprocessing_place)
                setting_str = setting_str + '\tProject data preprocessing result are stored in "' + \
                    data_preprocessing_place + '"\r'
                self.ui.textEdit.setText(setting_str)
            except:
                setting_str = setting_str + \
                    '\tCreating project data preprocessing storage folder ("' + \
                    fileplace + '") failed. Quit the analysis.' + '\r'
                self.ui.textEdit.setText(setting_str)
                return None
        currenttime = str(datetime.datetime.now()) + '\r'
        setting_str = setting_str + sharp_str+currenttime + 'Data preprocessing:' + '\r'
        self.ui.textEdit.setText(setting_str)
        QApplication.processEvents()
        # staring data preprocessing
        preprocessing_success = False
        control_label = self.ui.lineEdit_8.text()
        treat_label = self.ui.lineEdit_7.text()
        group_label = self.ui.lineEdit_5.text()
        permu_loops = self.ui.lineEdit_6.text()
        minimum_taxa_abundance_control = self.ui.lineEdit_9.text()

        try:
            qmd_dataPreprocessing_exe.preprocess_qmd(data_preprocessing_place, datafilename, minimum_taxa_detection_num, control_label, treat_label,
                                                     group_label, permu_loops, minimum_taxa_abundance_control, os_p)
            preprocessing_success = True
        except:
            pass

        taxaintomodel = []
        try:
            taxaintomodel = pd.read_csv(
                data_preprocessing_place+'/taxa_Into_Model.csv', header=0, index_col=0)
        except:
            preprocessing_success = False
        if preprocessing_success:
            if len(taxaintomodel) > 5:
                setting_str = setting_str + '\tData preprocessing success. ' + \
                    str(len(taxaintomodel))+' taxa into QMD analysis.' + '\r'
                self.ui.textEdit.setText(setting_str)
                QApplication.processEvents()
            if len(taxaintomodel) < 5:
                setting_str = setting_str + '\tAfter data preprocessing, ' + \
                    str(len(taxaintomodel)) + \
                    ' taxa left. The number is too small. Quit the analysis' + '\r'
                self.ui.textEdit.setText(setting_str)
                QApplication.processEvents()
                return None
        if not preprocessing_success:
            setting_str = setting_str + \
                '\tData preprocessing failed. Make sure the datafile follows the data template. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        currenttime = str(datetime.datetime.now()) + '\r'
        setting_str = setting_str + sharp_str+currenttime + 'QMD analysis:' + '\r'
        self.ui.textEdit.setText(setting_str)
        QApplication.processEvents()
        fdr = self.ui.comboBox_3.currentText()
        plot = self.ui.comboBox_5.currentText()
        confidence_level = [0.99, 0.95,
                            0.90][self.ui.comboBox_4.currentIndex()]
        qmd_optimization_exe.qmd_optimization(fileplace, data_preprocessing_place, control_label,
                                              treat_label, group_label, permu_loops, fdr, os_p, plot, confidence_level)
        analysis_sucess = False
        try:
            analysis_res = pd.read_csv(
                fileplace+'/QMD_result.csv', header=0, index_col=False)
            analysis_res_stacked = pd.read_csv(
                fileplace + '/QMD_result_stacked.csv', header=0, index_col=False)
            if len(analysis_res) > 0 and len(analysis_res_stacked) > 0:
                analysis_sucess = True
            currenttime = str(datetime.datetime.now()) + '\r'
            setting_str = setting_str + '\tThe quantified total microbial density changes between groups is ' + str(round(list(analysis_res.loc[analysis_res['type'] == 'qmd_total_density_change', 'value'])[0], 2)) + '.\r' +\
                '\tThe top positive changed taxa is '+str(analysis_res_stacked.loc[len(analysis_res_stacked)-1, 'item'])+' with '+str(round(analysis_res_stacked.loc[len(analysis_res_stacked)-1, 'qmd_density_diff'], 3)) + ' density fold increases.\r' + \
                '\tThe top negative changed taxa is ' + str(analysis_res_stacked.loc[0, 'item']) + ' with ' + str(round(analysis_res_stacked.loc[0, 'qmd_density_diff'], 3)) + ' density fold decreases.\r' + \
                '\tThe analysis result can be found in "'+fileplace+'"' + '\r'+ \
                sharp_str + currenttime + 'QMD analysis done.' + '\r'
            self.ui.textEdit.setText(setting_str)
            QApplication.processEvents()
            project_dict = {'projectname': self.ui.lineEdit_2.text(),
                            'permu': self.ui.lineEdit_6.text(),
                            'confidenceI': self.ui.comboBox_4.currentIndex(),
                            'fdr': self.ui.comboBox_3.currentIndex(),
                            'plot': self.ui.comboBox_5.currentIndex(),
                            'minimum_taxa_detection_samples': self.ui.lineEdit_4.text(),
                            'minimum_taxa_median_abundance': self.ui.lineEdit_9.text(),
                            'grouplabel': self.ui.lineEdit_5.text(),
                            'controllabel': self.ui.lineEdit_8.text(),
                            'treatlabel': self.ui.lineEdit_7.text(),
                            'datafile': self.ui.lineEdit_3.text(),
                            'loggs': setting_str
                            }
            load_save_project.save_dict(
                fileplace + '/project.qmd', project_dict)

            if os_p == 'Windows':
                logfileplace = (
                    fileplace + '/project.qmd').replace('/', '\\')
            else:
                logfileplace = (fileplace + '/projectname.qmd')
            setting_str = setting_str + 'The project configuration and logs are stored in "' + logfileplace + '."\r'
            self.ui.textEdit.setText(setting_str)
            QApplication.processEvents()
        except:
            pass
        if analysis_sucess == False:
            setting_str = setting_str + \
                'QMD analysis failed. Make sure the datafile follows the data template. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            QApplication.processEvents()
            return None


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = query_window()
    window.show()
    sys.exit(app.exec_())
