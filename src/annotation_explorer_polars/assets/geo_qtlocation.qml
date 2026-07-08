import QtQuick 2.15
import QtQuick.Controls 2.15
import QtLocation 6.0
import QtPositioning 6.0

Rectangle {
    id: root
    color: "#dbeafe"

    property var pointsModel: bridge ? bridge.points : []
    property var highlightedSamples: bridge ? bridge.highlightSamples : []
    property string statusText: bridge ? bridge.statusText : ""

    function isHighlighted(sampleName) {
        if (!highlightedSamples) {
            return false
        }
        return highlightedSamples.indexOf(String(sampleName)) !== -1
    }

    Plugin {
        id: osmPlugin
        name: "osm"
    }

    Map {
        id: mapView
        anchors.fill: parent
        plugin: osmPlugin
        center: QtPositioning.coordinate(18, 8)
        zoomLevel: 1.2
        copyrightsVisible: true

        Repeater {
            model: root.pointsModel

            delegate: MapQuickItem {
                coordinate: QtPositioning.coordinate(Number(modelData.lat), Number(modelData.lon))
                anchorPoint.x: marker.width / 2
                anchorPoint.y: marker.height / 2

                sourceItem: Rectangle {
                    id: marker
                    width: root.isHighlighted(modelData.sample) ? 18 : 12
                    height: width
                    radius: width / 2
                    color: root.isHighlighted(modelData.sample) ? "#f97316" : "#2563eb"
                    border.width: root.isHighlighted(modelData.sample) ? 3 : 1.6
                    border.color: root.isHighlighted(modelData.sample) ? "#7c2d12" : "#ffffff"
                    opacity: 0.92

                    ToolTip.visible: markerArea.containsMouse
                    ToolTip.text: {
                        var country = modelData.country ? String(modelData.country) + "\\n" : ""
                        return String(modelData.sample || "Sample") + "\\n" + country + "count: " + String(modelData.count || 0)
                    }

                    MouseArea {
                        id: markerArea
                        anchors.fill: parent
                        hoverEnabled: true
                    }
                }
            }
        }

        Timer {
            id: fitTimer
            interval: 120
            repeat: false
            onTriggered: {
                if (mapView.mapItems.length > 0) {
                    mapView.fitViewportToVisibleMapItems()
                    if (mapView.zoomLevel > 5.8) {
                        mapView.zoomLevel = 5.8
                    }
                } else {
                    mapView.center = QtPositioning.coordinate(18, 8)
                    mapView.zoomLevel = 1.2
                }
            }
        }

        Component.onCompleted: fitTimer.restart()
    }

    Rectangle {
        anchors.left: parent.left
        anchors.top: parent.top
        anchors.margins: 12
        width: Math.min(parent.width - 24, 420)
        radius: 10
        color: "#f8fafce8"
        border.color: "#cbd5e1"
        border.width: 1
        visible: true

        Column {
            anchors.fill: parent
            anchors.margins: 10
            spacing: 4

            Text {
                text: "Geo QML"
                color: "#0f172a"
                font.pixelSize: 13
                font.bold: true
            }

            Text {
                width: parent.width
                text: root.statusText
                wrapMode: Text.WordWrap
                color: "#334155"
                font.pixelSize: 12
            }
        }
    }

    Connections {
        target: bridge

        function onPointsChanged() {
            root.pointsModel = bridge.points
            fitTimer.restart()
        }

        function onHighlightSamplesChanged() {
            root.highlightedSamples = bridge.highlightSamples
        }

        function onStatusTextChanged() {
            root.statusText = bridge.statusText
        }
    }
}
