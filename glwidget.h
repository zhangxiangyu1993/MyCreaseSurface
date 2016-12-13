#ifndef GLWIDGET_H
#define GLWIDGET_H

#include "datadefine.h"
#include "DefineData.h"
#include <QGLWidget>
#include <QWheelEvent>
#include <QMouseEvent>
#include "GLUT.H"

class GLWidget : public QGLWidget
{
	Q_OBJECT
public:
	GLWidget();
	~GLWidget(void);

protected:
	void initializeGL();
	void resizeGL(int width, int height);
	void paintGL();
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent *event);
	void computeBound();
	void transformInputPoints();
	void drawAxis();
	void computeScale(GLfloat mv[4][4], GLfloat s[3], bool bColumn = true);
	void computeScale(GLfloat mv[16], GLfloat s[3], bool bColumn = true);
	void computeScale(GLdouble mv[16], GLdouble s[3], bool bColumn = true);
	void computeScale(GLdouble mv[4][4], GLdouble s[3], bool bColumn = true);
	void drawConstainedFaces();

private:	
	Vertices m_inputdata;
	Vertices m_origindata;
	Outputcells m_outputCells;

	GLfloat rotationX;
	GLfloat rotationY;
	GLfloat rotationZ;
	QPoint lastPos;

	GLfloat scrollAdd;
	GLfloat scroll;
	QFont*	    m_font;

	bool m_points, m_line;

	double m_minX,m_maxX,m_minY,m_maxY,m_minZ,m_maxZ;
};
#endif
