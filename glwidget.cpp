#include "glwidget.h"
#include <QDesktopWidget>
#include <QApplication>
#include <QFile>
#include <QTextStream>
#include <math.h>
#include <limits>
#include <assert.h>
#include <fstream>


using namespace std;

GLWidget::GLWidget()
{
	QDesktopWidget* desktop = QApplication::desktop();
	setGeometry(0, 0, 700, 700);

	move(100, 10);
	setWindowTitle( "My Cresase Surface" );
	//setFixedSize(700, 700);

	rotationX = 0.0;
	rotationY = 0.0;
	rotationZ = 0.0;
	scroll = 0.0;
	scrollAdd = 0.0;

	m_font = new QFont( "Times", 10);

	ifstream infilePoint("points.txt");
	ifstream infileOrder("pointorders.txt");
	ifstream infileOriginPoint("genvol.txt");

	/*int PointNum = 5032, TriNum = 2849;*/
	int PointNum = 0, TriNum = 0;
	int OriginNum = WSampleNum * HSampleNum * DSampleNum;
 	infilePoint >> PointNum;
	infileOrder >> TriNum;

	for (int i = 0; i < PointNum; i++)
	{
		Vertex t_vertex;
		infilePoint >> t_vertex.x >> t_vertex.y >> t_vertex.z;

		m_inputdata.push_back(t_vertex);
	}

	for (int j = 0; j < TriNum; j++)
	{
		Outputcell t_ocell;
		infileOrder >> t_ocell._v[0] >> t_ocell._v[1] >> t_ocell._v[2];

		m_outputCells.push_back(t_ocell);
	}

	for (int i = 0; i < OriginNum; i++)
	{
		Vertex t_vertex;

		infileOriginPoint >> t_vertex.proper;
		infileOriginPoint >> t_vertex.x >> t_vertex.y >> t_vertex.z;

		m_origindata.push_back(t_vertex);
	}

	m_points = true;
	m_line = true;

	m_minX = 0.0;
	m_maxX = 0.0;
	m_minY = 0.0;
	m_maxY = 0.0;
	m_minZ = 0.0;
	m_maxZ = 0.0;
// 	computeBound();
// 	transformInputPoints();
}

GLWidget::~GLWidget(void)
{
}

void GLWidget::initializeGL()
{
	glShadeModel(GL_SMOOTH);
	glClearColor(1.0, 1.0, 1.0, 0.5);
	glClearDepth(1.0);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
}

void GLWidget::resizeGL(int width, int height)
{
	if(height == 0)
		height = 1;
	glViewport(0, 0, (GLint)width, (GLint)height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//glOrtho(m_minX,m_maxX,m_minY,m_maxY,m_minZ,m_minZ);
	gluPerspective(45.0, (GLfloat)width/(GLfloat)height, 0.1, 100.0);
	glMatrixMode(GL_MODELVIEW);

}

void GLWidget::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
#if 1
	glPushMatrix();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, -4.0f + scroll);
	glRotatef(rotationX, 1.0, 0.0, 0.0);
	glRotatef(rotationY, 0.0, 1.0, 0.0);
	glRotatef(rotationZ, 0.0, 0.0, 1.0);
	glEnable(GL_POINT_SMOOTH);

	if (m_points)
	{
		glColor3f(1.0f, 0.0f, 0.0f);
		glPointSize(2);
		for(size_t i = 0;i < m_outputCells.size(); ++i)
		{
			int _v[4];
			for (int j = 0;j < 4;++j)
			{
				_v[j] = m_outputCells[i]._v[j];
			}

			glBegin(GL_POINTS);
			glVertex3f(m_inputdata[_v[0]].x, m_inputdata[_v[0]].y, m_inputdata[_v[0]].z);
			glVertex3f(m_inputdata[_v[1]].x, m_inputdata[_v[1]].y, m_inputdata[_v[1]].z);
			glVertex3f(m_inputdata[_v[2]].x, m_inputdata[_v[2]].y, m_inputdata[_v[2]].z);
			glEnd();
	
		}

	
		glPointSize(2);

		glBegin(GL_POINTS);
		for (size_t i = 0; i < m_origindata.size();++i)
		{
			if (m_origindata[i].proper > 0.0)
			{
				glColor3f(0.0f, 0.0f, 1.0f);
				glVertex3f(m_origindata[i].x, m_origindata[i].y, m_origindata[i].z);
			}
// 			}else{
// 				glColor3f(0.0f, 0.0f, 0.0f);
// 			}

			
		}

		glEnd();

	}

	if (m_line)
	{
		glColor3f(0.0f, 0.0f, 0.0f);
		glLineWidth(1);
	
		for(size_t i = 0;i < m_outputCells.size();++i)
		{
			int _v[4];
			for (int j = 0;j < 4;++j)
			{
				_v[j] = m_outputCells[i]._v[j];
			}
			glBegin(GL_LINE_LOOP);
			glVertex3f(m_inputdata[_v[0]].x, m_inputdata[_v[0]].y, m_inputdata[_v[0]].z);
			glVertex3f(m_inputdata[_v[1]].x, m_inputdata[_v[1]].y, m_inputdata[_v[1]].z);
			glVertex3f(m_inputdata[_v[2]].x, m_inputdata[_v[2]].y, m_inputdata[_v[2]].z);
			glEnd();
			
		}

	}

	glPopMatrix();
#endif
	drawAxis();
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
	lastPos = event->pos();
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
	GLfloat dx = GLfloat(event->x() - lastPos.x()) / width();
	GLfloat dy = GLfloat(event->y() - lastPos.y()) / height();
	if(event->buttons() & Qt::LeftButton)
	{
		rotationX += 180 * dy;
		rotationZ += 180 * dx;
		updateGL();
	}
	lastPos = event->pos();
}

void GLWidget::wheelEvent(QWheelEvent *event)
{
	double numDegrees = event->delta() / 8.0;
	scrollAdd = numDegrees * 0.03;
	scroll += scrollAdd;
	event->accept();
	updateGL();
}

void GLWidget::computeBound()
{
	m_minX = m_inputdata[0].x;
	m_maxX = m_inputdata[0].x;
	m_minY = m_inputdata[0].y;
	m_maxY = m_inputdata[0].y;
	m_minZ = m_inputdata[0].z;
	m_maxZ = m_inputdata[0].z;

	for (int i = 1;i < m_inputdata.size();++i)
	{
		if ( m_inputdata[i].x < m_minX ) m_minX = m_inputdata[i].x;
		if ( m_inputdata[i].x > m_maxX ) m_maxX = m_inputdata[i].x;
		if ( m_inputdata[i].y < m_minY ) m_minY = m_inputdata[i].y;
		if ( m_inputdata[i].y > m_maxY ) m_maxY = m_inputdata[i].y;
		if ( m_inputdata[i].z < m_minZ ) m_minZ = m_inputdata[i].z;
		if ( m_inputdata[i].z > m_maxZ ) m_maxZ = m_inputdata[i].z;
	}

}

void GLWidget::transformInputPoints(){
	if (m_inputdata.size() <= 0)
	{
		return;
	}

	Vertices _transformInputdata;
	//初始化设置输入点不需要进行坐标转换
	bool pointTransformOrNot = false;
	for (int i = 0;i < m_inputdata.size();++i)
	{
		if ((abs(m_inputdata[i].x)>1) || (abs(m_inputdata[i].y)>1) || (abs(m_inputdata[i].z)>1))
		{
			pointTransformOrNot = true;
			break;
		}
	}	
	if (pointTransformOrNot)
	{
		//转换
		_transformInputdata.resize(m_inputdata.size());
		for (int i = 0;i < m_inputdata.size();++i)
		{
			Vertex _ver;
			_ver.x = (m_inputdata[i].x - (m_maxX+m_minX)/2.0) / ((m_maxX-m_minX)/2.0);
			_ver.y = (m_inputdata[i].y - (m_maxY+m_minY)/2.0) / ((m_maxY-m_minY)/2.0);
			_ver.z = (m_inputdata[i].z - (m_maxZ+m_minZ)/2.0) / ((m_maxZ-m_minZ)/2.0);
			_ver.proper = m_inputdata[i].proper;
			_transformInputdata[i] = _ver;
		}
		
		m_inputdata = _transformInputdata;
	}

	_transformInputdata.resize(0);

	pointTransformOrNot = false;
	for (int i = 0;i < m_origindata.size();++i)
	{
		if ((abs(m_origindata[i].x)>1) || (abs(m_origindata[i].y)>1) || (abs(m_origindata[i].z)>1))
		{
			pointTransformOrNot = true;
			break;
		}
	}

	if (pointTransformOrNot)
	{
		//转换
		_transformInputdata.resize(m_origindata.size());
		for (int i = 0;i < m_origindata.size();++i)
		{
			Vertex _ver;
			_ver.x = (m_origindata[i].x - (m_maxX+m_minX)/2.0) / ((m_maxX-m_minX)/2.0);
			_ver.y = (m_origindata[i].y - (m_maxY+m_minY)/2.0) / ((m_maxY-m_minY)/2.0);
			_ver.z = (m_origindata[i].z - (m_maxZ+m_minZ)/2.0) / ((m_maxZ-m_minZ)/2.0);
			_ver.proper = m_origindata[i].proper;
			_transformInputdata[i] = _ver;
		}

		m_origindata = _transformInputdata;
	}
	
}

void GLWidget::drawAxis(void)
{
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glViewport(0, 0, 100, 100);	 
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(-1.5, 1.5, -1.5, 1.5, 0.1, 10);
	glTranslatef(0, 0, -1.1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	GLdouble mv[4][4], s[3];
	glGetDoublev(GL_MODELVIEW_MATRIX, (GLdouble *)mv);
	mv[3][0] = mv[3][1] = mv[3][2] = 0.0;
	mv[3][3] = 1.0;
	glLoadMatrixd((GLdouble *)mv);
	computeScale(mv, s);
	glScalef(1.0 / s[0], 1.0 / s[1], 1.0 / s[2]);

	glRotatef(rotationX, 1.0, 0.0, 0.0);
	glRotatef(rotationY, 0.0, 1.0, 0.0);
	glRotatef(rotationZ, 0.0, 0.0, 1.0);

	glColor3f(0.0,0.0,0.0);
	renderText(1.1f, 0, 0, "X",*m_font,2000);
	renderText(0, 1.1f, 0, "Y",*m_font,2000);
	renderText(0, 0, 1.1f, "Z",*m_font,2000);

	glColor3f(0.0, 0.0, 0.0);
	glutSolidSphere(0.1,12,12);  
	//绘制小坐标
	glColor3f(1.0, 0.0, 0.0);
	glLineWidth(3);
	glBegin(GL_LINES);
	{
		glColor3f(1.0, 0.0, 0.0);
		glVertex3f(0.0f, 0, 0); // x axis
		glVertex3f(1.0f, 0, 0);

		glVertex3f(1.0f, 0, 0);
		glVertex3f(0.8f, 0.2f, 0);
		glVertex3f(1.0f, 0, 0);
		glVertex3f(0.8f, -.2f, 0);

		glColor3f(0.0, 1.0, 0.0);
		glVertex3f(0, 0.0f, 0); // y axis
		glVertex3f(0, 1.0f, 0);

		glVertex3f(0, 1.0f, 0);
		glVertex3f(0.2f, 0.8f, 0);
		glVertex3f(0, 1.0f, 0);
		glVertex3f(-0.2f, 0.8f, 0);

		glColor3f(0.0, 0.0, 1.0);
		glVertex3f(0, 0, 0.0f); // z axis
		glVertex3f(0, 0, 1.0f);

		glVertex3f(0, 0, 1.0f);
		glVertex3f(0, 0.2f, 0.8f);
		glVertex3f(0, 0, 1.0f);
		glVertex3f(0,-0.2f, 0.8f);
	}
	glEnd();
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glPopAttrib();
}

/*****************************************************
功能描述:按比例缩放
参数：
返回：
作者：
******************************************************/
void GLWidget::computeScale(GLfloat mv[4][4], GLfloat s[3], bool bColumn)
{
	GLdouble alpha[4][4], beta[3];

	for(size_t i = 0; i < 4; i++)
	{
		for(size_t j = 0; j < 4; j ++)
			alpha[i][j] = mv[i][j];
	}

	computeScale(alpha, beta);
	copy(beta, beta + 3, s);
}

void GLWidget::computeScale(GLfloat mv[16], GLfloat s[3], bool bColumn)
{
	GLdouble alpha[4][4], beta[3];

	for(size_t i = 0; i < 4; i++)
	{
		for(size_t j = 0; j < 4; j ++)
			alpha[i][j] = mv[i*4 + j];
	}
	computeScale(alpha, beta); 
	copy(beta, beta + 3, s);
}

void GLWidget::computeScale(GLdouble mv[16], GLdouble s[3], bool bColumn)
{
	GLdouble alpha[4][4];
	for(size_t i = 0; i < 4; i++)
	{
		for(size_t j = 0; j < 4; j ++)
			alpha[i][j] = mv[i*4 + j];
	}
	computeScale(alpha, s);
}

void GLWidget::computeScale(GLdouble mv[4][4], GLdouble s[3], bool bColumn /*= true*/)
{
	for(int i = 0; i < 3; ++i)
	{
		s[i] =0;
		for(int j = 0; j < 3; ++j)
		{
			if(bColumn)
			{
				s[i] += mv[i][j] * mv[i][j];
			}
			else
			{
				s[i] += mv[j][i] * mv[j][i];
			}
		}
		assert(s[i] >= 0);
		if(fabs(s[i]) <= numeric_limits<GLfloat>::epsilon())
		{
			s[i] = numeric_limits<float>::epsilon();
		}
		s[i] = sqrt(s[i]);
	}
}